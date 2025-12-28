import pytest
import math
from satellite_visibility.core import (
    j2000_to_cirs_precession_only,
    geodetic_to_itrs,
    itrs_to_cirs,
    calculate_elevation,
    check_visibility
)


class TestJ2000ToCirsPrecession:
    """Тесты для функции прецессии J2000 → CIRS"""
    
    def test_identity_at_j2000_epoch(self):
        """В эпоху J2000 (t=0) преобразование должно быть близко к единичной матрице"""
        pos_j2000 = [6378137, 0, 0]
        pos_cirs = j2000_to_cirs_precession_only(pos_j2000, 0)
        diff = sum((pos_cirs[i] - pos_j2000[i])**2 for i in range(3))**0.5
        assert diff < 1000, f"Разность {diff} м должна быть <1000 м"
    
    def test_vector_length_preservation(self):
        """Проверка сохранения длины вектора (ортогональность матрицы)"""
        pos_j2000 = [7000000, 1000000, 2000000]
        pos_cirs = j2000_to_cirs_precession_only(pos_j2000, 7305)  # 20 лет
        
        len_j2000 = sum(x**2 for x in pos_j2000)**0.5
        len_cirs = sum(x**2 for x in pos_cirs)**0.5
        
        assert abs(len_j2000 - len_cirs) < 0.1, f"Разность длин {abs(len_j2000 - len_cirs)} м должна быть <0.1 м"
    
    def test_precession_magnitude_100_years(self):
        """Проверка величины прецессии за 100 лет"""
        pos_j2000 = [6378137, 0, 0]
        pos_cirs_100y = j2000_to_cirs_precession_only(pos_j2000, 36525)  # 100 лет
        
        angle = math.atan2(pos_cirs_100y[1], pos_cirs_100y[0])
        angle_arcsec = angle * 180 * 3600 / math.pi
        
        assert 4000 < angle_arcsec < 6000, f"Прецессия {angle_arcsec}\" должна быть ~5000\""


class TestGeodeticToItrs:
    """Тесты для конвертации геодезических координат в ITRS"""
    
    def test_equator_point(self):
        """Точка на экваторе должна дать правильный радиус"""
        pos = geodetic_to_itrs(0, 0, 0)
        expected_radius = 6378137.0
        actual_radius = (pos[0]**2 + pos[1]**2 + pos[2]**2)**0.5
        
        assert abs(actual_radius - expected_radius) < 1, f"Радиус {actual_radius} м должен быть {expected_radius} м"
    
    def test_north_pole(self):
        """Северный полюс должен дать правильные координаты"""
        pos_pole = geodetic_to_itrs(90, 0, 0)
        
        assert abs(pos_pole[0]) < 1, f"X координата {pos_pole[0]} должна быть ≈ 0"
        assert abs(pos_pole[1]) < 1, f"Y координата {pos_pole[1]} должна быть ≈ 0"
        assert 6356000 < pos_pole[2] < 6357000, f"Z координата {pos_pole[2]} должна быть ≈ 6356752 м"
    
    def test_height_difference(self):
        """Проверка разности высот"""
        pos_sea = geodetic_to_itrs(45, 0, 0)
        pos_1km = geodetic_to_itrs(45, 0, 1000)
        
        height_diff = ((pos_1km[0]**2 + pos_1km[1]**2 + pos_1km[2]**2)**0.5 - 
                       (pos_sea[0]**2 + pos_sea[1]**2 + pos_sea[2]**2)**0.5)
        
        assert 990 < height_diff < 1010, f"Разность высот {height_diff} м должна быть ~1000 м"


class TestItrsToCircsRotation:
    """Тесты для поворота ITRS → CIRS"""
    
    def test_era_daily_change(self):
        """ERA за день должна изменяться на ~361°"""
        era_0 = 2 * math.pi * (0.7790572732640 + 1.00273781191135448 * 0)
        era_1 = 2 * math.pi * (0.7790572732640 + 1.00273781191135448 * 1)
        
        era_diff = era_1 - era_0
        era_diff_deg = math.degrees(era_diff)
        
        assert 360 < era_diff_deg < 362, f"Изменение ERA {era_diff_deg}° должно быть ~361°"
    
    def test_coordinate_rotation(self):
        """Поворот координат за день"""
        pos_itrs = [6378137, 0, 0]
        pos_cirs_0 = itrs_to_cirs(pos_itrs, 0)
        pos_cirs_1 = itrs_to_cirs(pos_itrs, 1)
        
        func_angle_0 = math.atan2(pos_cirs_0[1], pos_cirs_0[0])
        func_angle_1 = math.atan2(pos_cirs_1[1], pos_cirs_1[0])
        func_diff = math.degrees(func_angle_1 - func_angle_0)
        
        assert 360 < abs(func_diff) < 362, f"Поворот {func_diff}° должен быть ~361°"


class TestSatelliteVisibility:
    """Тесты видимости спутника"""
    
    def test_satellite_at_zenith(self):
        """Спутник в зените над северным полюсом должен быть виден"""
        sat_pos_j2000 = [0, 0, 7000000]
        elevation, visible = check_visibility(sat_pos_j2000, 90, 0, 0, 0)
        
        assert elevation > 80, f"Угол возвышения {elevation}° должен быть >80°"
        assert visible, "Спутник должен быть виден"
    
    def test_elevation_threshold(self):
        """Проверка порога видимости 15 градусов"""
        sat_pos_j2000 = [0, 0, 7000000]
        elevation, visible = check_visibility(sat_pos_j2000, 90, 0, 0, 0, 15)
        
        assert visible, "Спутник должен быть виден при пороге 15°"
    
    def test_satellite_behind_horizon(self):
        """Спутник на противоположной стороне Земли не должен быть виден"""
        sat_pos_j2000 = [-7000000, 0, 0]
        elevation, visible = check_visibility(sat_pos_j2000, 0, 0, 0, 0)
        
        assert elevation < 0, f"Угол возвышения {elevation}° должен быть отрицательным"
        assert not visible, "Спутник не должен быть виден"


class TestCalculateElevation:
    """Тесты вычисления угла возвышения"""
    
    def test_elevation_calculation_zenith(self):
        """Тест вычисления угла для спутника в зените"""
        sat_pos_cirs = [0, 0, 7000000]
        obs_pos_cirs = [0, 0, 6378137]
        obs_lat_rad = math.radians(90)
        
        elevation = calculate_elevation(sat_pos_cirs, obs_pos_cirs, obs_lat_rad)
        
        assert elevation > 80, f"Угол возвышения {elevation}° должен быть близок к 90°"
    
    def test_elevation_returns_number(self):
        """Функция должна возвращать числовое значение"""
        sat_pos_cirs = [1000000, 0, 0]
        obs_pos_cirs = [6378137, 0, 0]
        obs_lat_rad = 0
        
        elevation = calculate_elevation(sat_pos_cirs, obs_pos_cirs, obs_lat_rad)
        
        assert isinstance(elevation, (int, float)), "Результат должен быть числом"