import numpy as np
import matplotlib.pyplot as plt
from .core import (
    j2000_to_cirs_precession_only,
    geodetic_to_itrs,
    itrs_to_cirs,
    calculate_elevation
)


def visualize_satellite_visibility(sat_pos_j2000, obs_lat_deg, obs_lon_deg, obs_height_m, 
                                 jd_days_from_j2000, min_elevation_deg=0):
    """
    Визуализация видимости спутника
    
    Args:
        sat_pos_j2000: [x, y, z] координаты спутника в J2000 (метры)
        obs_lat_deg: широта наблюдателя в градусах (float)
        obs_lon_deg: долгота наблюдателя в градусах (float)
        obs_height_m: высота наблюдателя в метрах (float)
        jd_days_from_j2000: дни от эпохи J2000 (float)
        min_elevation_deg: минимальный угол возвышения в градусах (float)
    
    Returns:
        (elevation_deg, is_visible): угол возвышения и видимость (float, bool)
    """
    # Вычисляем позиции
    sat_pos_cirs = j2000_to_cirs_precession_only(sat_pos_j2000, jd_days_from_j2000)
    obs_pos_itrs = geodetic_to_itrs(obs_lat_deg, obs_lon_deg, obs_height_m)
    obs_pos_cirs = itrs_to_cirs(obs_pos_itrs, jd_days_from_j2000)
    
    elevation_deg = calculate_elevation(sat_pos_cirs, obs_pos_cirs, np.radians(obs_lat_deg))
    is_visible = elevation_deg > min_elevation_deg
    
    # Создаем фигуру
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Находим плоскость через центр Земли, наблюдателя и спутник
    earth_center = np.array([0, 0, 0])
    observer = np.array(obs_pos_cirs)
    satellite = np.array(sat_pos_cirs)
    
    # Проецируем на плоскость для 2D визуализации
    # Используем наблюдателя как базовую ось
    obs_dist = np.linalg.norm(observer)
    sat_dist = np.linalg.norm(satellite)
    
    # Угол между наблюдателем и спутником от центра Земли
    cos_angle = np.dot(observer, satellite) / (obs_dist * sat_dist)
    angle = np.arccos(np.clip(cos_angle, -1, 1))
    
    # 2D координаты в плоскости
    obs_2d = np.array([obs_dist, 0])
    sat_2d = np.array([sat_dist * np.cos(angle), sat_dist * np.sin(angle)])
    earth_2d = np.array([0, 0])
    
    # Радиус Земли
    earth_radius = 6378137
    
    # Рисуем Землю
    circle = plt.Circle(earth_2d, earth_radius, color='lightblue', alpha=0.7, label='Земля')
    ax.add_patch(circle)
    
    # Рисуем точки
    ax.plot(*earth_2d, 'ko', markersize=8, label='Центр Земли')
    ax.plot(*obs_2d, 'ro', markersize=8, label='Наблюдатель')
    ax.plot(*sat_2d, 'bs', markersize=8, label='Спутник')
    
    # Рисуем линии
    ax.plot([earth_2d[0], obs_2d[0]], [earth_2d[1], obs_2d[1]], 'r--', alpha=0.5)
    ax.plot([obs_2d[0], sat_2d[0]], [obs_2d[1], sat_2d[1]], 'g-', linewidth=2, label='Линия визирования')
    
    # Зона видимости (касательные к Земле от наблюдателя)
    if obs_dist > earth_radius:
        # Угол к горизонту
        horizon_angle = np.arcsin(earth_radius / obs_dist)
        
        # Касательные точки
        tangent_angle1 = -horizon_angle
        tangent_angle2 = horizon_angle
        
        # Длина касательной
        tangent_length = np.sqrt(obs_dist**2 - earth_radius**2)
        
        # Точки касания
        tang1_x = obs_2d[0] + tangent_length * np.cos(tangent_angle1)
        tang1_y = obs_2d[1] + tangent_length * np.sin(tangent_angle1)
        
        tang2_x = obs_2d[0] + tangent_length * np.cos(tangent_angle2)
        tang2_y = obs_2d[1] + tangent_length * np.sin(tangent_angle2)
        
        # Рисуем зону видимости
        ax.plot([obs_2d[0], tang1_x], [obs_2d[1], tang1_y], 'orange', alpha=0.7, label='Граница видимости')
        ax.plot([obs_2d[0], tang2_x], [obs_2d[1], tang2_y], 'orange', alpha=0.7)
        
        # Заливаем зону видимости
        max_dist = max(sat_dist, obs_dist) * 1.2
        
        # Создаем веер лучей в зоне видимости
        angles = np.linspace(tangent_angle1, tangent_angle2, 20)
        for angle in angles[::3]:  # Каждый третий луч
            end_x = obs_2d[0] + max_dist * np.cos(angle)
            end_y = obs_2d[1] + max_dist * np.sin(angle)
            ax.plot([obs_2d[0], end_x], [obs_2d[1], end_y], 'yellow', alpha=0.3, linewidth=0.5)
    
    # Настройки графика
    ax.set_xlim(-earth_radius * 0.5, max(sat_dist, obs_dist) * 1.2)
    ax.set_ylim(-earth_radius * 1.5, earth_radius * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Заголовок с результатами
    visibility_text = "ВИДЕН" if is_visible else "НЕ ВИДЕН"
    color = 'green' if is_visible else 'red'
    
    ax.set_title(f'Видимость спутника\n'
                f'Угол возвышения: {elevation_deg:.1f}°\n'
                f'Статус: {visibility_text}', 
                fontsize=14, color=color, weight='bold')
    
    ax.set_xlabel('Расстояние (м)')
    ax.set_ylabel('Расстояние (м)')
    
    plt.tight_layout()
    plt.show()
    
    return elevation_deg, is_visible


def check_visibility_with_plot(sat_pos_j2000, obs_lat_deg, obs_lon_deg, obs_height_m, 
                              jd_days_from_j2000, min_elevation_deg=0):
    """
    Проверка видимости спутника с визуализацией
    
    Args:
        sat_pos_j2000: [x, y, z] координаты спутника в J2000 (метры)
        obs_lat_deg: широта наблюдателя в градусах (float)
        obs_lon_deg: долгота наблюдателя в градусах (float)
        obs_height_m: высота наблюдателя в метрах (float)
        jd_days_from_j2000: дни от эпохи J2000 (float)
        min_elevation_deg: минимальный угол возвышения в градусах (float)
    
    Returns:
        (elevation_deg, is_visible): угол возвышения и видимость (float, bool)
    """
    elevation_deg, is_visible = visualize_satellite_visibility(
        sat_pos_j2000, obs_lat_deg, obs_lon_deg, obs_height_m, 
        jd_days_from_j2000, min_elevation_deg
    )
    
    visibility_text = "виден" if is_visible else "не виден"
    print(f"Угол возвышения: {elevation_deg:.1f} градусов, КА {visibility_text}")
    
    return elevation_deg, is_visible