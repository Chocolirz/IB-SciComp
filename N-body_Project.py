import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Body:
    """表示一个具有质量、位置和速度的物体"""
    
    def __init__(self, mass, position, velocity, color=None, name=None):
        """
        初始化一个物体
        
        参数:
            mass: 质量
            position: 位置 (2D numpy数组)
            velocity: 速度 (2D numpy数组)
            color: 可视化颜色
            name: 物体名称
        """
        self.mass = mass
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.color = color if color else np.random.rand(3,)
        self.name = name if name else f"Body_{id(self)}"
        self.trajectory = [self.position.copy()]
        
    def update_velocity(self, acceleration, dt):
        """根据加速度更新速度"""
        self.velocity += acceleration * dt
        
    def update_position(self, dt):
        """根据速度更新位置"""
        self.position += self.velocity * dt
        self.trajectory.append(self.position.copy())
        
    def distance_to(self, other_body):
        """计算到另一个物体的距离"""
        return np.linalg.norm(self.position - other_body.position)
        
    def __repr__(self):
        return f"{self.name}: m={self.mass}, pos={self.position}, vel={self.velocity}"


class GravitationalSystem:
    """表示一个引力系统，包含多个物体"""
    
    def __init__(self, G=6.67430e-11, softening=0.1):
        """
        初始化引力系统
        
        参数:
            G: 万有引力常数
            softening: 软化参数，防止距离过近时引力过大
        """
        self.G = G
        self.softening = softening
        self.bodies = []
        self.time = 0
        
    def add_body(self, body):
        """向系统中添加一个物体"""
        self.bodies.append(body)
        
    def calculate_acceleration(self, body):
        """计算指定物体受到的总加速度"""
        acceleration = np.zeros(2)
        
        for other in self.bodies:
            if body is not other:
                # 计算相对位置向量
                r_vector = other.position - body.position
                # 计算距离（加入软化参数防止除零错误）
                distance = np.linalg.norm(r_vector) + self.softening
                # 计算引力加速度
                acceleration += self.G * other.mass * r_vector / (distance ** 3)
                
        return acceleration
    
    def update_system(self, dt):
        """更新整个系统一个时间步长"""
        # 首先计算所有物体的加速度
        accelerations = [self.calculate_acceleration(body) for body in self.bodies]
        
        # 然后更新所有物体的速度和位置
        for body, acceleration in zip(self.bodies, accelerations):
            body.update_velocity(acceleration, dt)
            body.update_position(dt)
            
        self.time += dt
        
    def simulate(self, total_time, dt, progress_callback=None):
        """模拟系统一段时间"""
        num_steps = int(total_time / dt)
        
        for step in range(num_steps):
            self.update_system(dt)
            
            if progress_callback and step % 10 == 0:
                progress_callback(step / num_steps)
                
    def get_kinetic_energy(self):
        """计算系统的总动能"""
        return sum(0.5 * body.mass * np.linalg.norm(body.velocity) ** 2 for body in self.bodies)
    
    def get_potential_energy(self):
        """计算系统的总势能"""
        potential_energy = 0
        n = len(self.bodies)
        
        for i in range(n):
            for j in range(i + 1, n):
                distance = self.bodies[i].distance_to(self.bodies[j])
                potential_energy -= self.G * self.bodies[i].mass * self.bodies[j].mass / distance
                
        return potential_energy
    
    def get_total_energy(self):
        """计算系统的总能量"""
        return self.get_kinetic_energy() + self.get_potential_energy()


class SimulationVisualizer:
    """可视化引力系统模拟"""
    
    def __init__(self, system, trail_length=100):
        self.system = system
        self.trail_length = trail_length
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.ax.set_xlabel('X Position')
        self.ax.set_ylabel('Y Position')
        self.ax.set_title('N-Body Gravitational Simulation')
        self.ax.grid(True)
        
        # 初始化散点图和轨迹线
        self.scatter = self.ax.scatter([], [], s=50)
        self.trails = []
        
        # 为每个物体创建轨迹线
        for body in system.bodies:
            trail, = self.ax.plot([], [], '-', color=body.color, alpha=0.3)
            self.trails.append(trail)
            
        # 设置坐标轴范围
        self.ax.set_xlim(-15, 15)
        self.ax.set_ylim(-15, 15)
        self.ax.set_aspect('equal', 'box')
        
    def init_animation(self):
        """初始化动画"""
        self.scatter.set_offsets(np.empty((0, 2)))
        for trail in self.trails:
            trail.set_data([], [])
        return [self.scatter] + self.trails
    
    def update_animation(self, frame):
        """更新动画帧"""
        # 更新系统
        self.system.update_system(0.01)
        
        # 更新散点图
        positions = np.array([body.position for body in self.system.bodies])
        self.scatter.set_offsets(positions)
        
        # 更新轨迹
        for i, body in enumerate(self.system.bodies):
            trail_x = [pos[0] for pos in body.trajectory[-self.trail_length:]]
            trail_y = [pos[1] for pos in body.trajectory[-self.trail_length:]]
            self.trails[i].set_data(trail_x, trail_y)
        '''
        # 动态调整坐标轴范围
        all_positions = np.concatenate([body.trajectory[-self.trail_length:] for body in self.system.bodies])
        if len(all_positions) > 0:
            x_min, y_min = np.min(all_positions, axis=0)
            x_max, y_max = np.max(all_positions, axis=0)
            margin = 0.1
            self.ax.set_xlim(x_min - margin, x_max + margin)
            self.ax.set_ylim(y_min - margin, y_max + margin)
        '''
        return [self.scatter] + self.trails
    
    def animate(self, dt=0.01, frames=1000, interval=20):
        """创建并运行动画"""
        self.animation = FuncAnimation(
            self.fig, self.update_animation, frames=frames,
            init_func=self.init_animation, blit=True, interval=interval
        )
        plt.show()


# 示例用法
if __name__ == "__main__":
    # 创建引力系统
    sun_mass = 1000
    jupiter_mass = 1
    jupiter_avg_distance = 10
    system = GravitationalSystem(G=1.0)  # 归一化G值以简化计算
    G = system.G  # 使用系统的G值
    orbital_velocity = np.sqrt(G * sun_mass / jupiter_avg_distance)
    
    # 创建一些示例物体，这里用的是太阳、木星和一个小行星，初始位于木星领先相位的拉格朗日点上。
    # 但是，这里我们假设了太阳初始是不动的，会导致模型不甚准确，不过，小星星仍然会围绕拉格朗日点做运动，
    # 只不过运动比较复杂，详见刘川老师所写的《理论力学》一书。
    planet1 = Body(mass=sun_mass, position=[0, 0], velocity=[0, 0], color='red', name="SUN")
    planet2 = Body(mass=jupiter_mass, position=[0, jupiter_avg_distance], velocity=[orbital_velocity, 0], color='green', name="JUPITER")
    planet3 = Body(mass=jupiter_mass/1000, position=[0.5*np.sqrt(3)*jupiter_avg_distance, 0.5*jupiter_avg_distance], velocity=[0.5*orbital_velocity, -0.5*np.sqrt(3)*orbital_velocity], color='blue', name="Planet")
    
    # 将物体添加到系统
    system.add_body(planet1)
    system.add_body(planet2)
    system.add_body(planet3)
    
    # 创建可视化器并运行动画
    visualizer = SimulationVisualizer(system, trail_length=200)
    visualizer.animate(dt=0.01, frames=1000, interval=20)
    
    # 如果需要批量模拟而不是动画，可以使用以下代码：
    # system.simulate(total_time=10, dt=0.01)
    # 然后可以访问每个物体的轨迹进行分析