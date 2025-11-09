#include <cmath>
#include <cstring>
#include "structures.h"
#include "util.h"

const double Pi = 3.14159265358979323846;

// ============================================
// 初始化单圈传输矩阵
// ============================================
void calc_transfer_matrix(double tunex, double tuney, double* Mtransfer)
{
    // 检查指针
    if (Mtransfer == nullptr) {
        error_exit("Mtransfer is nullptr in calc_transfer_matrix", __FILE__, __LINE__);
    }
    
    // 初始化为单位矩阵
    memset(Mtransfer, 0, 36 * sizeof(double));
    for (int i = 0; i < 6; ++i) {
        Mtransfer[i * 6 + i] = 1.0;
    }
    
    // 横向 x 方向的旋转矩阵（前2×2块）
    Mtransfer[0 * 6 + 0] = cos(2.0 * Pi * tunex);
    Mtransfer[0 * 6 + 1] = sin(2.0 * Pi * tunex);
    Mtransfer[1 * 6 + 0] = -sin(2.0 * Pi * tunex);
    Mtransfer[1 * 6 + 1] = cos(2.0 * Pi * tunex);
    
    // 横向 y 方向的旋转矩阵（中间2×2块）
    Mtransfer[2 * 6 + 2] = cos(2.0 * Pi * tuney);
    Mtransfer[2 * 6 + 3] = sin(2.0 * Pi * tuney);
    Mtransfer[3 * 6 + 2] = -sin(2.0 * Pi * tuney);
    Mtransfer[3 * 6 + 3] = cos(2.0 * Pi * tuney);
    
    // 纵向部分保持为单位矩阵（z, delta）
    // Mtransfer[4*6+4] = 1.0;  // 已经在初始化时设置
    // Mtransfer[5*6+5] = 1.0;
}

// ============================================
// 初始化束流
// ============================================
void initialize_beam(Beam* beam, const BeamParams* params)
{
    // 复制基本参数
    beam->nbunch = params->nbunch;
    beam->nmp = params->nmp;
    beam->nmp0 = params->nmp;
    beam->macrosize = params->macrosize;
    beam->energy = params->energy;
    
    beam->sigma_x = params->sigma_x;
    beam->sigma_y = params->sigma_y;
    beam->sigma_z = params->sigma_z;
    
    beam->tunex = params->tunex;
    beam->tuney = params->tuney;
    beam->tunes = params->tunes;
    
    beam->tracking_mode = params->tracking_mode;
    beam->U0 = params->U0;
    
    // 计算相对论因子
    double electron_mass_eV = 0.5109989461e6;  // 电子静止能量 (eV)
    beam->gamma = beam->energy / electron_mass_eV;
    beam->beta = sqrt(1.0 - 1.0 / (beam->gamma * beam->gamma));
    
    // 分配粒子内存
    allocate_particle_memory(beam);
    
    // 生成粒子分布
    generate_particle_distribution(beam);
    
    // 计算束流参数
    compute_beam_parameters(beam);
    
    // 初始化传输矩阵（如果需要横向追踪）
    if (beam->tracking_mode == 1) {
        calc_transfer_matrix(beam->tunex, beam->tuney, beam->Mtransfer);
    }
}

// ============================================
// 单圈传输映射
// ============================================
void apply_one_turn_map(Beam* beam, const SimulationParams* params)
{
    double tRev = params->CIRC / (beam->beta * get_light_speed());
    double eta = params->eta;
    double c = get_light_speed();
    
    int total_particles = beam->nbunch * beam->nmp;
    
    if (beam->tracking_mode == 0) {
        // ============ 模式 0: 仅纵向追踪 ============
        apply_longitudinal_map_only(beam, tRev, eta, c);
    }
    else if (beam->tracking_mode == 1) {
        // ============ 模式 1: 横向+纵向追踪 ============
        apply_full_6d_map(beam, tRev, eta, c);
    }
    else {
        error_exit("Unknown tracking_mode in apply_one_turn_map", __FILE__, __LINE__);
    }
}

// ============================================
// 仅纵向传输映射
// ============================================
void apply_longitudinal_map_only(Beam* beam, double tRev, double eta, double c)
{
    int total_particles = beam->nbunch * beam->nmp;
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // 仅更新纵向坐标
        // 注意：这里的 coord[4] 可能是 z 或 t，根据实际定义
        // 假设 coord[4] = z（纵向位置）, coord[5] = delta（能量偏差）
        
        double delta = coord[5];
        
        // 纵向传输公式：
        // Δt = tRev × (1 + delta × η)
        // Δz = c × Δt = c × tRev × (1 + delta × η)
        coord[4] += c * tRev * (1.0 + delta * eta);
        
        // 横向坐标保持不变（tracking_mode = 0）
        // coord[0], coord[1], coord[2], coord[3] 不变
    }
}

// ============================================
// 横向+纵向传输映射（6D）
// ============================================
void apply_full_6d_map(Beam* beam, double tRev, double eta, double c)
{
    int total_particles = beam->nbunch * beam->nmp;
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // ===== 横向传输：使用旋转矩阵 =====
        double x_old = coord[0];
        double px_old = coord[1];
        double y_old = coord[2];
        double py_old = coord[3];
        
        // 新坐标 = Mtransfer × 旧坐标
        coord[0] = beam->Mtransfer[0 * 6 + 0] * x_old + 
                   beam->Mtransfer[0 * 6 + 1] * px_old;
        coord[1] = beam->Mtransfer[1 * 6 + 0] * x_old + 
                   beam->Mtransfer[1 * 6 + 1] * px_old;
        
        coord[2] = beam->Mtransfer[2 * 6 + 2] * y_old + 
                   beam->Mtransfer[2 * 6 + 3] * py_old;
        coord[3] = beam->Mtransfer[3 * 6 + 2] * y_old + 
                   beam->Mtransfer[3 * 6 + 3] * py_old;
        
        // ===== 纵向传输：使用色散关系 =====
        double delta = coord[5];
        coord[4] += c * tRev * (1.0 + delta * eta);
    }
}

// ============================================
// 更新单个粒子的纵向坐标（独立函数）
// ============================================
void update_longitudinal_coordinate(double* particle_coord, 
                                   double tRev, double eta, double c)
{
    // particle_coord[4] = z
    // particle_coord[5] = delta
    
    double delta = particle_coord[5];
    particle_coord[4] += c * tRev * (1.0 + delta * eta);
}

// ============================================
// 同步辐射能损
// ============================================
void apply_synchrotron_radiation(Beam* beam)
{
    int total_particles = beam->nbunch * beam->nmp;
    double U0_normalized = beam->U0 / beam->energy;  // 归一化能损
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // delta 减少
        coord[5] -= U0_normalized;
    }
}

// ============================================
// 辐射阻尼（可选）
// ============================================
void apply_radiation_damping(Beam* beam, double damping_time)
{
    // 简化模型：指数阻尼
    // damping_time 单位：圈数
    
    int total_particles = beam->nbunch * beam->nmp;
    double damping_factor = exp(-1.0 / damping_time);
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // 阻尼作用于所有坐标（简化）
        for (int dim = 0; dim < 6; ++dim) {
            coord[dim] *= damping_factor;
        }
    }
}

// ============================================
// 量子激发（可选）
// ============================================
void apply_quantum_excitation(Beam* beam, double sigma_delta)
{
    // 量子激发：能量方向的随机kick
    
    int total_particles = beam->nbunch * beam->nmp;
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // 仅作用于 delta（能量偏差）
        coord[5] += gaussian_random(0.0, sigma_delta);
    }
}

// ============================================
// 粒子管理工具函数
// ============================================
double get_particle_coordinate(const Beam* beam, int particle_id, int dim)
{
    if (particle_id < 0 || particle_id >= beam->nbunch * beam->nmp) {
        error_exit("Invalid particle_id in get_particle_coordinate", 
                   __FILE__, __LINE__);
    }
    if (dim < 0 || dim >= DIM_COORD) {
        error_exit("Invalid dimension in get_particle_coordinate", 
                   __FILE__, __LINE__);
    }
    
    return beam->coords[particle_id * DIM_COORD + dim];
}

void set_particle_coordinate(Beam* beam, int particle_id, int dim, double value)
{
    if (particle_id < 0 || particle_id >= beam->nbunch * beam->nmp) {
        error_exit("Invalid particle_id in set_particle_coordinate", 
                   __FILE__, __LINE__);
    }
    if (dim < 0 || dim >= DIM_COORD) {
        error_exit("Invalid dimension in set_particle_coordinate", 
                   __FILE__, __LINE__);
    }
    
    beam->coords[particle_id * DIM_COORD + dim] = value;
}

int get_particle_bunch_id(int particle_id, int nmp)
{
    return particle_id / nmp;
}

int determine_bucket(double z, double CIRC, int h)
{
    // 将 z 坐标映射到 bucket ID (0 到 h-1)
    // 假设 z 的范围是 [0, CIRC)
    
    double bucket_length = CIRC / h;
    
    // 归一化到 [0, CIRC)
    double z_normalized = fmod(z, CIRC);
    if (z_normalized < 0) z_normalized += CIRC;
    
    int bucket_id = (int)(z_normalized / bucket_length);
    
    // 边界处理
    if (bucket_id >= h) bucket_id = h - 1;
    if (bucket_id < 0) bucket_id = 0;
    
    return bucket_id;
}

// ============================================
// 计算束流参数
// ============================================
void compute_beam_parameters(Beam* beam)
{
    // 每束团粒子数
    beam->np = beam->macrosize * beam->nmp;
    
    // 平均束流（假设均匀分布在h个bucket中）
    // IbDC = nbunch × np × e × fRev
    // 这里需要从其他参数计算 fRev，暂时留空
    // 实际应该在 Wake.cpp 或 Cavity.cpp 中计算
    
    // beam->IbDC = ...;  // 需要 fRev
}

// ============================================
// 分配粒子内存
// ============================================
void allocate_particle_memory(Beam* beam)
{
    int total_size = beam->nbunch * beam->nmp * DIM_COORD;
    beam->coords = new double[total_size];
    
    // 初始化为0
    memset(beam->coords, 0, total_size * sizeof(double));
}

// ============================================
// 释放粒子内存
// ============================================
void free_particle_memory(Beam* beam)
{
    if (beam->coords != nullptr) {
        delete[] beam->coords;
        beam->coords = nullptr;
    }
}

// ============================================
// 生成粒子分布（高斯分布）
// ============================================
void generate_particle_distribution(Beam* beam)
{
    int total_particles = beam->nbunch * beam->nmp;
    
    for (int i = 0; i < total_particles; ++i) {
        double* coord = &beam->coords[i * DIM_COORD];
        
        // 横向分布（高斯）
        coord[0] = gaussian_random(0.0, beam->sigma_x);  // x
        coord[1] = gaussian_random(0.0, beam->sigma_x);  // px（简化：用sigma_x）
        coord[2] = gaussian_random(0.0, beam->sigma_y);  // y
        coord[3] = gaussian_random(0.0, beam->sigma_y);  // py
        
        // 纵向分布（高斯）
        coord[4] = gaussian_random(0.0, beam->sigma_z);  // z
        coord[5] = gaussian_random(0.0, beam->sigma_z / 100.0);  // delta（简化）
        
        // 其他维度初始化
        coord[6] = 1.0;  // liveflag = 1（存活）
        coord[7] = -1.0; // loss_iele = -1（未损失）
        coord[8] = -1.0; // loss_iturn = -1
        coord[9] = (double)i;  // init_index
    }
}

// ============================================
// 计算束流统计量
// ============================================
void compute_beam_statistics(const Beam* beam, BeamStats* stats)
{
    int total_particles = beam->nbunch * beam->nmp;
    
    // 初始化
    for (int dim = 0; dim < 6; ++dim) {
        stats->mean[dim] = 0.0;
        stats->std[dim] = 0.0;
    }
    
    // 计算均值
    for (int i = 0; i < total_particles; ++i) {
        const double* coord = &beam->coords[i * DIM_COORD];
        for (int dim = 0; dim < 6; ++dim) {
            stats->mean[dim] += coord[dim];
        }
    }
    for (int dim = 0; dim < 6; ++dim) {
        stats->mean[dim] /= total_particles;
    }
    
    // 计算标准差
    for (int i = 0; i < total_particles; ++i) {
        const double* coord = &beam->coords[i * DIM_COORD];
        for (int dim = 0; dim < 6; ++dim) {
            double diff = coord[dim] - stats->mean[dim];
            stats->std[dim] += diff * diff;
        }
    }
    for (int dim = 0; dim < 6; ++dim) {
        stats->std[dim] = sqrt(stats->std[dim] / total_particles);
    }
}
