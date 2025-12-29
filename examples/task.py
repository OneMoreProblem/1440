from satellite_visibility import *

position_eci = [4435144, -2137297, 4670064]
observation_time = 8084.185608609847
observer_lat = 45.920266
observer_lon = -64.342286
min_elevation = 15.0

elevation_deg, visibility = check_visibility(position_eci, observer_lat, observer_lon, 0, 
                                            observation_time, min_elevation_deg=0)
visibility = 'виден' if visibility else 'не виден'

print('Угол возвышения:', elevation_deg, 'градусов, КА', visibility)