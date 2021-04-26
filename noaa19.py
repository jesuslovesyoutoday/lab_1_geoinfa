import sgp4
from math import radians
from math import degrees
import math
import numpy as np
from datetime import datetime, timedelta
import datetime as dt
import matplotlib.pyplot as plt
from pyorbital.orbital import Orbital
import spacetrack
import requests
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time


"""_______________________________Parcing_TLE_______________________________"""

def from_strings(tle_file, sat_name):
    r = requests.get(tle_file, allow_redirects=True)
    open('tle.txt', 'wb').write(r.content)
    txt = open('tle.txt', 'r').read().split('\n')
    for i in range(len(txt)):
        if (txt[i] == sat_name):
            print('#_____________________________TLE_DATA:______________________________#')
            print(sat_name, txt[i + 1], txt[i + 2], sep='\n')
            print(' ')
            tle = [sat_name, txt[i + 1], txt[i + 2]]
            return tle

"""____________________________ECI_Coordinates______________________________"""

def eci_coords():

	#orb elements from TLE
    dolg_vosh_uzl = 247.4627
    arg_per = 130.5360
    ist_anom = 325.0288
    nakl_orb = 51.6416
    sredn_dvij = 15.72125391
    grav_potential = 398603*10e9
    eccentricity = 0.0006703
    u = arg_per + ist_anom

    #big poluos
    a = (((1/sredn_dvij)**2)*grav_potential)**(1/3)
    
    #focalniy parametr
    p = (1 - (eccentricity)**2)*a
    
    #radius-vector
    r = p/(1 + eccentricity* math.cos(ist_anom*math.pi/180))

    x = r*(math.cos(dolg_vosh_uzl*math.pi/180)*math.cos(u*math.pi/180) - math.sin(dolg_vosh_uzl*math.pi/180)*math.sin(u*math.pi/180)*math.cos(nakl_orb*math.pi/180))

    y = r*(math.sin(dolg_vosh_uzl*math.pi/180)*math.cos(u*math.pi/180) + math.cos(dolg_vosh_uzl*math.pi/180)*math.sin(u*math.pi/180)*math.cos(nakl_orb*math.pi/180))

    z = r*math.sin(u*math.pi/180)*math.sin(nakl_orb*math.pi/180)

    return [x, y, z]
    
"""_______________________ECEF_Coordinates_from_ECI__________________________"""

def ecef_coords(xyz, now):

    cartrep = coord.CartesianRepresentation(*xyz, unit=u.m)
    gcrs = coord.GCRS(cartrep, obstime=now) 
    itrs = gcrs.transform_to(coord.ITRS(obstime=now)) 
    loc = coord.EarthLocation(*itrs.cartesian.xyz) 
    print('') 
    print(loc.lat, loc.lon, loc.height)
    return [loc.lat, loc.lon, loc.height]
    
"""_____________________Input_required_time_period___________________________"""
    
def define_period():

	start_date_entry = input('Enter a start date (YYYY-MM-DD-HH-MM)  :   ')
	start_year, start_month, start_day, start_hour, start_minute = map(int, start_date_entry.split('-'))
	start_date = datetime(start_year, start_month, start_day, start_hour, start_minute)

	finish_date_entry = input('Enter a finish date (YYYY-MM-DD-HH-MM) :   ')
	finish_year, finish_month, finish_day, finish_hour, finish_minute = map(int, finish_date_entry.split('-'))
	finish_date = datetime(finish_year, finish_month, finish_day, finish_hour, finish_minute)
	
	print(' ')

	utc_start_date = start_date - timedelta(hours=3)
	utc_finish_date = finish_date - timedelta(hours=3)
	
	return [utc_start_date, utc_finish_date];
	
"""_____________________________LK_coordinates______________________________"""

def get_LK_coords():

	LK_lat = 55.93013
	LK_lon = 37.51832
	height_over_sea_level = 0.2
	Earth_radius = 6378.137
	distance_to_LK = Earth_radius + height_over_sea_level
	LK_x = distance_to_LK * math.cos(radians(LK_lat)) * math.cos(radians(LK_lon))
	LK_y = distance_to_LK * math.sin(radians(LK_lon)) * math.cos(radians(LK_lat))
	LK_z = distance_to_LK * math.sin(radians(LK_lat))
	return [LK_x, LK_y, LK_z]
	
"""_______________________Get_azimuths_and_elevations_______________________"""

def SHIT_WHAT(LK_coords):

	time = define_period();
	sat_x = []
	sat_y = []
	sat_z = []
	
	all_thetas = []
	all_phis = []
	all_times = []
	
	Earth_radius = 6378.137
    
	while (time[0] != time[1]):
    
		orb = Orbital("N", line1=tle[1], line2=tle[2])
		s_deg_lon, s_deg_lat, s_alt = orb.get_lonlatalt(time[0])
		s_rad_lon = radians(s_deg_lon)
		s_rad_lat = radians(s_deg_lat)
		s_r = s_alt + Earth_radius
		s_x = s_r * math.cos(s_rad_lat) * math.cos(s_rad_lon)
		s_y = s_r * math.sin(s_rad_lon) * math.cos(s_rad_lat)
		s_z = s_r * math.sin(s_rad_lat)
    	
		sat_x.append(s_x)
		sat_y.append(s_y)
		sat_z.append(s_z)
    	
		LK = LK_coords[0] ** 2 + LK_coords[1] ** 2 + LK_coords[2] ** 2
		sk = LK_coords[0]* s_x + LK_coords[1] * s_y + LK_coords[2] * s_z
    
		dist_to_fl = (sk - LK)/(LK ** 0.5)
		dist_to_sat = ((s_x - LK_coords[0]) ** 2 + (s_y - LK_coords[1]) ** 2 + (s_z - LK_coords[2]) ** 2) ** 0.5
		theta = math.asin(dist_to_fl / dist_to_sat)
		theta = degrees(theta)
		phi = 0
    	
		if theta >= 0:
			x_n = 0
			y_n = 0
			z_n = LK / LK_coords[2]
			vec_nor = [x_n - LK_coords[0], y_n - LK_coords[1], z_n - LK_coords[2]]
			len_nor = (vec_nor[0] ** 2 + vec_nor[1] ** 2 + vec_nor[2] ** 2) ** 0.5
			t = (LK - sk)/LK
			x_p = s_x + t * LK_coords[0]
			y_p = s_y + t * LK_coords[1]
			z_p = s_z + t * LK_coords[2]
			vec_pr = [x_p - LK_coords[0], y_p - LK_coords[1], z_p - LK_coords[2]]
			len_pr = (vec_pr[0] ** 2 + vec_pr[1] ** 2 + vec_pr[2] ** 2) ** 0.5
	
			phi = math.acos((vec_nor[0] * vec_pr[0] + vec_nor[1] * vec_pr[1] + vec_nor[2] * vec_pr[2]) / (len_nor * len_pr))
			phi = degrees(phi)
	
			x1 = vec_nor[0]
			y1 = vec_nor[1]
			z1 = vec_nor[2]
			vec_est = [y1 * LK_coords[2] - z1 * LK_coords[1],
	                   -(x1 * LK_coords[2] - LK_coords[0] * z1),
	                   x1 * LK_coords[1] - y1 * LK_coords[0]]
			len_est = (vec_est[0] ** 2 + vec_est[1] ** 2 + vec_est[2] ** 2) ** 0.5
			phi1 = math.acos((vec_est[0] * vec_pr[0] + vec_est[1] * vec_pr[1] + vec_est[2] * vec_pr[2]) / (len_est * len_pr))
			phi1 = degrees(phi1)

			if phi1 > 90:
				phi = -phi + 360
		
		all_thetas.append(theta)
		all_phis.append(phi)
		all_times.append(time[0] + timedelta(hours = 3))
		time[0] += timedelta(minutes = 1)
   
	events_thetas = []
	events_phis = []
	events_times = []

	cur_event_thetas = []
	cur_event_phis = []
	cur_event_times = []
    
	for i in range(len(all_thetas)):
		if all_thetas[i] >= 0:
			cur_event_thetas.append(all_thetas[i])
			cur_event_phis.append(radians(all_phis[i]))
			cur_event_times.append(all_times[i])
		else:
			if cur_event_thetas:
				events_thetas.append(cur_event_thetas)
				events_phis.append(cur_event_phis)
				events_times.append(cur_event_times)
			cur_event_thetas = []
			cur_event_phis = []
			cur_event_times = []
	return [events_times, events_phis, events_thetas, sat_x, sat_y, sat_z]       		
    
"""_____________________________Final_graphs________________________________"""
    
def graph():
	
	LK_coords = get_LK_coords();
	mega_list = SHIT_WHAT(LK_coords);
	
	for i in range(len(mega_list[0])):
	
		print('Time period:  ', mega_list[0][i][0], '  --  ', mega_list[0][i][len(mega_list[0][i]) - 1])
		print('Azimuth: ', str(degrees(mega_list[1][i][0]))[:-12], 'degrees; ', end=' ')
		print('Max elevation: ', str(max(mega_list[2][i]))[:-12], 'degrees')

	sf = plt.figure()
	ax = sf.add_subplot(111, projection='3d')
	ax.plot(mega_list[3], mega_list[4], mega_list[5])
	ax.scatter(LK_coords[0], LK_coords[1], LK_coords[2], color='red')
	sf.set_size_inches(7, 7)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	ax.set_theta_zero_location('N')
	ax.set_theta_direction(-1)
	ax.set_rlim(bottom=90, top=0)
	for phi, theta in zip(mega_list[1], mega_list[2]):
		ax.plot(phi, theta)
	fig.set_size_inches(7, 7)
	plt.show()    	
    
    
"""__________________________________MAIN___________________________________"""    

    
tle_file = "https://celestrak.com/NORAD/elements/active.txt"

satellite = "NOAA 19                 "

tle = from_strings(tle_file, satellite)

graph()
