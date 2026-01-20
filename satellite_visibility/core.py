import math
import numpy as np


def j2000_to_cirs_precession_only(pos_j2000, jd_days_from_j2000):
    """
    Конвертация ECI J2000 в CIRS с учетом только прецессии (IAU 2006)
    
    Args:
        pos_j2000: [x, y, z] координаты в системе J2000 (метры)
        jd_days_from_j2000: дни от эпохи J2000 (float)
    
    Returns:
        [x, y, z] координаты в системе CIRS (метры)
    """
    T = jd_days_from_j2000 / 36525.0

    # Углы в миллисекундах дуги по IAU 2006 (прецессия без нутации)
    psi_mas = (5038.481507 * T +
               0.00000014 * T**2 -
               0.000000013 * T**3)
    omega_mas = (84381.406 -
                 46.836769 * T -
                 0.00001831 * T**2 +
                 0.00200340 * T**3 -
                 5.76e-7 * T**4)
    chi_mas = (10.556403 * T -
               0.000010 * T**2 -
               3.281e-7 * T**3 +
               0.00000014 * T**4)

    # В радианы
    arcsec_to_rad = math.pi / (180 * 3600)
    psi = psi_mas * arcsec_to_rad
    omega = omega_mas * arcsec_to_rad
    chi = chi_mas * arcsec_to_rad

    # Матрицы вращения
    def R1(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])

    def R3(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    # Матрица прецессии IAU 2006 (без нутации)
    P = R3(chi) @ R1(-omega) @ R3(-psi) @ R1(omega)

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

def calculate_elevation(sat_pos_cirs, obs_pos_cirs):
    """
    Вычисление угла возвышения спутника над горизонтом
    
    Args:
        sat_pos_cirs: [x, y, z] позиция спутника в CIRS (метры)
        obs_pos_cirs: [x, y, z] позиция наблюдателя в CIRS (метры)
    
    Returns:
        угол возвышения в градусах (float)
    """
    # Вектор от наблюдателя к спутнику
    sat_obs_vec = np.array(sat_pos_cirs) - np.array(obs_pos_cirs)
    
    # Локальная вертикаль (направление из центра Земли к наблюдателю)
    z_local = np.array(obs_pos_cirs)
    z_local /= np.linalg.norm(z_local)  # Нормируем
    
    # Проекция вектора на вертикаль
    dot_product = np.dot(sat_obs_vec, z_local)
    vec_length = np.linalg.norm(sat_obs_vec)
    
    # Угол возвышения
    elevation_rad = math.asin(dot_product / vec_length)
    return math.degrees(elevation_rad)

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
    
    elevation_deg = calculate_elevation(
        sat_pos_cirs, obs_pos_cirs
    )
    
    is_visible = elevation_deg > min_elevation_deg
    
    return elevation_deg, is_visible