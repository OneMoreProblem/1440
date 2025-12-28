import math
import numpy as np

def j2000_to_cirs_precession_only(pos_j2000, jd_days_from_j2000):
    """
    Конвертация ECI J2000 в CIRS с учетом только прецессии
    
    Args:
        pos_j2000: [x, y, z] координаты в системе J2000 (метры)
        jd_days_from_j2000: дни от эпохи J2000 (float)
    
    Returns:
        [x, y, z] координаты в системе CIRS (метры)
    """
    T = jd_days_from_j2000 / 36525.0
    
    psi = (-0.041775 + 5038.481484 * T + 1.5584175 * T**2) * math.pi / (180 * 3600)
    phi = (84381.412819 - 46.811016 * T + 0.0511268 * T**2) * math.pi / (180 * 3600)
    gamma = (-0.052928 + 10.556378 * T + 0.4932044 * T**2) * math.pi / (180 * 3600)
    eps0 = 84381.406 * math.pi / (180 * 3600)
    
    def R1(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])
    
    def R3(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    
    P = R3(-gamma) @ R1(phi) @ R3(-psi) @ R1(-eps0)
    return (P @ np.array(pos_j2000)).tolist()

def geodetic_to_itrs(lat_deg, lon_deg, height_m):
    """
    Конвертация геодезических координат в ITRS (WGS84)
    
    Args:
        lat_deg: широта в градусах (float)
        lon_deg: долгота в градусах (float)
        height_m: высота над эллипсоидом в метрах (float)
    
    Returns:
        [x, y, z] координаты в системе ITRS (метры)
    """
    a = 6378137.0
    e2 = 0.00669438002290
    
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    
    N = a / math.sqrt(1 - e2 * math.sin(lat_rad)**2)
    
    x = (N + height_m) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (N + height_m) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (N * (1 - e2) + height_m) * math.sin(lat_rad)
    
    return [x, y, z]

def itrs_to_cirs(pos_itrs, jd_days_from_j2000):
    """
    Конвертация ITRS в CIRS с учетом вращения Земли
    
    Args:
        pos_itrs: [x, y, z] координаты в системе ITRS (метры)
        jd_days_from_j2000: дни от эпохи J2000 (float)
    
    Returns:
        [x, y, z] координаты в системе CIRS (метры)
    """
    era = 2 * math.pi * (0.7790572732640 + 1.00273781191135448 * jd_days_from_j2000)
    
    c = math.cos(-era)
    s = math.sin(-era)
    
    x_cirs = c * pos_itrs[0] + s * pos_itrs[1]
    y_cirs = -s * pos_itrs[0] + c * pos_itrs[1]
    z_cirs = pos_itrs[2]
    
    return [x_cirs, y_cirs, z_cirs]

def calculate_elevation(sat_pos_cirs, obs_pos_cirs, obs_lat_rad):
    """
    Вычисление угла возвышения спутника над горизонтом
    
    Args:
        sat_pos_cirs: [x, y, z] позиция спутника в CIRS (метры)
        obs_pos_cirs: [x, y, z] позиция наблюдателя в CIRS (метры)
        obs_lat_rad: широта наблюдателя в радианах (float)
    
    Returns:
        угол возвышения в градусах (float)
    """
    sat_obs_vec = np.array(sat_pos_cirs) - np.array(obs_pos_cirs)
    vec_length = np.linalg.norm(sat_obs_vec)
    
    obs_pos = np.array(obs_pos_cirs)
    z_local = obs_pos / np.linalg.norm(obs_pos)
    
    dot_product = np.dot(sat_obs_vec, z_local)
    elevation_rad = math.asin(dot_product / vec_length)
    elevation_deg = math.degrees(elevation_rad)
    
    return elevation_deg

def check_visibility(sat_pos_j2000, obs_lat_deg, obs_lon_deg, obs_height_m, 
                    jd_days_from_j2000, min_elevation_deg=0):
    """
    Проверка видимости спутника
    
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
    sat_pos_cirs = j2000_to_cirs_precession_only(sat_pos_j2000, jd_days_from_j2000)
    obs_pos_itrs = geodetic_to_itrs(obs_lat_deg, obs_lon_deg, obs_height_m)
    obs_pos_cirs = itrs_to_cirs(obs_pos_itrs, jd_days_from_j2000)
    
    elevation_deg = calculate_elevation(sat_pos_cirs, obs_pos_cirs, math.radians(obs_lat_deg))
    is_visible = elevation_deg > min_elevation_deg
    
    return elevation_deg, is_visible