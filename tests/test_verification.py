import pytest
import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    CartesianRepresentation,
    ICRS,
    CIRS,
    AltAz,
    EarthLocation,
    GCRS
)
from astropy import units as u

# Импортируем самописные функции
from satellite_visibility.core import (
    j2000_to_cirs_precession_only,
    geodetic_to_itrs,
    itrs_to_cirs,
    calculate_elevation,
    check_visibility
)


class TestAstropyVerification:
    """Тесты для верификации самописных функций через сравнение с Astropy"""

    @pytest.fixture
    def test_data(self):
        """Фикстура с тестовыми данными"""
        return {
            'sat_pos_j2000': [4435144, -2137297, 4670064],
            'jd_days': 8084.185608609847,
            'obs_lat': 45.920266,
            'obs_lon': -64.342286,
            'obs_h': 0
        }

    def astropy_full_solution(self, sat_pos_j2000, obs_lat, obs_lon, obs_h, jd_days):
        """Решение с использованием Astropy (эталонное)"""
        time = Time("J2000") + jd_days * u.day
        
        # Создаем геоцентрические координаты спутника в GCRS (геоцентрическая система)
        cart = CartesianRepresentation(sat_pos_j2000 * u.m)
        gcrs = GCRS(cart, obstime=time)
        
        # Преобразуем в CIRS
        cirs_sat = gcrs.transform_to(CIRS(obstime=time))
        
        # Создаем местоположение наблюдателя
        location = EarthLocation.from_geodetic(obs_lon * u.deg, obs_lat * u.deg, obs_h * u.m)
        
        # Преобразуем положение наблюдателя в CIRS
        obs_itrs = location.get_itrs(time)
        obs_cirs = obs_itrs.transform_to(CIRS(obstime=time))
        
        # Вычисляем вектор от наблюдателя к спутнику в CIRS
        sat_obs_vector = CartesianRepresentation(
            cirs_sat.cartesian.xyz - obs_cirs.cartesian.xyz
        )
        
        # Создаем AltAz систему координат относительно наблюдателя
        altaz_frame = AltAz(obstime=time, location=location)
        
        # Создаем координаты для вектора от наблюдателя к спутнику
        topo = CIRS(sat_obs_vector, obstime=time, location=location).transform_to(altaz_frame)
        
        return topo.alt.degree

    def self_written_full_solution(self, sat_pos_j2000, obs_lat, obs_lon, obs_h, jd_days):
        """Решение с использованием самописных функций"""
        sat_pos_cirs = j2000_to_cirs_precession_only(sat_pos_j2000, jd_days)
        obs_pos_itrs = geodetic_to_itrs(obs_lat, obs_lon, obs_h)
        obs_pos_cirs = itrs_to_cirs(obs_pos_itrs, jd_days)
        
        elevation_deg = calculate_elevation(sat_pos_cirs, obs_pos_cirs)
        return elevation_deg

    def test_elevation_comparison(self, test_data):
        """Сравнение углов возвышения между Astropy и самописными функциями"""
        elev_astropy = self.astropy_full_solution(
            test_data['sat_pos_j2000'],
            test_data['obs_lat'],
            test_data['obs_lon'],
            test_data['obs_h'],
            test_data['jd_days']
        )
        
        elev_self = self.self_written_full_solution(
            test_data['sat_pos_j2000'],
            test_data['obs_lat'],
            test_data['obs_lon'],
            test_data['obs_h'],
            test_data['jd_days']
        )
        
        # Допустимая погрешность 0.1 градуса
        assert abs(elev_astropy - elev_self) < 0.1, \
            f"Разница углов возвышения {abs(elev_astropy - elev_self):.6f}° превышает допустимую погрешность 0.1°"

    def test_j2000_to_cirs_precession_only(self, test_data):
        """Сравнение преобразования J2000 → CIRS"""
        # Astropy решение
        time = Time("J2000") + test_data['jd_days'] * u.day
        cart = CartesianRepresentation(test_data['sat_pos_j2000'] * u.m)
        gcrs = GCRS(cart, obstime=time)
        cirs_astropy = gcrs.transform_to(CIRS(obstime=time))
        cirs_astropy_array = cirs_astropy.cartesian.xyz.to(u.m).value
        
        # Самописное решение
        cirs_self = np.array(j2000_to_cirs_precession_only(
            test_data['sat_pos_j2000'], 
            test_data['jd_days']
        ))
        
        # Сравнение с допустимой погрешностью 100 метров (влияние прецессии без нутации)
        diff = np.linalg.norm(cirs_astropy_array - cirs_self)
        assert diff < 100, \
            f"Разница в преобразовании J2000 → CIRS составляет {diff:.2f} м, что превышает допустимую погрешность 100 м"

    def test_geodetic_to_itrs(self, test_data):
        """Сравнение преобразования геодезических координат в ITRS"""
        # Astropy решение
        location = EarthLocation.from_geodetic(
            test_data['obs_lon'] * u.deg, 
            test_data['obs_lat'] * u.deg, 
            test_data['obs_h'] * u.m
        )
        obs_itrs_astropy = location.get_itrs()
        obs_itrs_astropy_array = obs_itrs_astropy.cartesian.xyz.to(u.m).value
        
        # Самописное решение
        obs_itrs_self = np.array(geodetic_to_itrs(
            test_data['obs_lat'],
            test_data['obs_lon'], 
            test_data['obs_h']
        ))
        
        # Сравнение с допустимой погрешностью 1 метр
        diff = np.linalg.norm(obs_itrs_astropy_array - obs_itrs_self)
        assert diff < 1, \
            f"Разница в преобразовании geodetic → ITRS составляет {diff:.2f} м, что превышает допустимую погрешность 1 м"

    def test_itrs_to_cirs(self, test_data):
        """Сравнение преобразования ITRS → CIRS"""
        # Получаем позицию наблюдателя в ITRS
        obs_pos_itrs = geodetic_to_itrs(
            test_data['obs_lat'],
            test_data['obs_lon'], 
            test_data['obs_h']
        )
        
        # Astropy решение
        time = Time("J2000") + test_data['jd_days'] * u.day
        location = EarthLocation.from_geodetic(
            test_data['obs_lon'] * u.deg, 
            test_data['obs_lat'] * u.deg, 
            test_data['obs_h'] * u.m
        )
        obs_itrs = location.get_itrs(time)
        obs_cirs_astropy = obs_itrs.transform_to(CIRS(obstime=time))
        obs_cirs_astropy_array = obs_cirs_astropy.cartesian.xyz.to(u.m).value
        
        # Самописное решение
        obs_cirs_self = np.array(itrs_to_cirs(obs_pos_itrs, test_data['jd_days']))
        
        # Сравнение с допустимой погрешностью 1 метр
        diff = np.linalg.norm(obs_cirs_astropy_array - obs_cirs_self)
        assert diff < 1, \
            f"Разница в преобразовании ITRS → CIRS составляет {diff:.2f} м, что превышает допустимую погрешность 1 м"

    def test_multiple_scenarios(self):
        """Тестирование нескольких сценариев для повышения надежности"""
        test_scenarios = [
            {
                'sat_pos_j2000': [4435144, -2137297, 4670064],
                'jd_days': 8084.185608609847,
                'obs_lat': 45.920266,
                'obs_lon': -64.342286,
                'obs_h': 0
            },
            {
                'sat_pos_j2000': [6000000, 2000000, 1000000],
                'jd_days': 10000.5,
                'obs_lat': 0,
                'obs_lon': 0,
                'obs_h': 1000
            },
            {
                'sat_pos_j2000': [-5000000, 3000000, -2000000],
                'jd_days': 15000.25,
                'obs_lat': -30,
                'obs_lon': 45,
                'obs_h': 500
            }
        ]
        
        for i, scenario in enumerate(test_scenarios):
            elev_astropy = self.astropy_full_solution(
                scenario['sat_pos_j2000'],
                scenario['obs_lat'],
                scenario['obs_lon'],
                scenario['obs_h'],
                scenario['jd_days']
            )
            
            elev_self = self.self_written_full_solution(
                scenario['sat_pos_j2000'],
                scenario['obs_lat'],
                scenario['obs_lon'],
                scenario['obs_h'],
                scenario['jd_days']
            )
            
            # Допустимая погрешность 0.1 градуса
            assert abs(elev_astropy - elev_self) < 0.1, \
                f"Сценарий {i+1}: Разница углов возвышения {abs(elev_astropy - elev_self):.6f}° превышает допустимую погрешность 0.1°"

    def test_visibility_check_consistency(self, test_data):
        """Проверка согласованности функции check_visibility с отдельными вычислениями"""
        # Используем отдельные функции
        sat_pos_cirs = j2000_to_cirs_precession_only(
            test_data['sat_pos_j2000'], 
            test_data['jd_days']
        )
        obs_pos_itrs = geodetic_to_itrs(
            test_data['obs_lat'],
            test_data['obs_lon'], 
            test_data['obs_h']
        )
        obs_pos_cirs = itrs_to_cirs(obs_pos_itrs, test_data['jd_days'])
        elevation_direct = calculate_elevation(sat_pos_cirs, obs_pos_cirs)
        
        # Используем функцию высшего уровня
        elevation_function, _ = check_visibility(
            test_data['sat_pos_j2000'],
            test_data['obs_lat'],
            test_data['obs_lon'],
            test_data['obs_h'],
            test_data['jd_days'],
            0
        )
        
        # Результаты должны совпадать
        assert abs(elevation_direct - elevation_function) < 1e-10, \
            f"Результаты check_visibility и отдельных функций не совпадают: {abs(elevation_direct - elevation_function)}"