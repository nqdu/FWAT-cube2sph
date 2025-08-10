import glob 

# Input and output files
input_files = glob.glob("STATIONS_*")
output_file = 'output.txt'
print(input_files)

# Read and combine all lines
lines = []
for fname in input_files:
    with open(fname) as f:
        lines.extend(f.readlines())

# Sort and deduplicate (lexically)
unique_sorted = sorted(set(line.strip() for line in lines))

print(unique_sorted)

# # Write to output.txt
# with open(output_file, 'w') as f:
#     for line in unique_sorted:
#         f.write(line + '\n')


# import numpy as np 
# def _geod2geoc(geographic_lat_deg, flattening=0.0033528106647474805):
#     """
#     Convert geographic (geodetic) latitude to geocentric latitude.
    
#     Parameters:
#         geographic_lat_deg : float
#             Geographic latitude in degrees.
#         flattening : float
#             Flattening of the ellipsoid (default is WGS-84).
    
#     Returns:
#         geocentric_lat_deg : float
#             Geocentric latitude in degrees.
#     """
#     phi = np.radians(geographic_lat_deg)
#     geocentric_phi = np.arctan((1 - flattening) ** 2 * np.tan(phi))
#     return np.degrees(geocentric_phi)

# def cal_dist_az_baz(lat1:float,lon1:float,lat2:float,lon2:float):
#     from obspy.geodetics import calc_vincenty_inverse

#     dist_in_m,az,baz = calc_vincenty_inverse(_geod2geoc(lat1),lon1,_geod2geoc(lat2),lon2,6371000.,0.)
    
#     return dist_in_m,az,baz

# print(cal_dist_az_baz(30,40,60,50))
# print(cal_dist_az_baz(60,50,30,40))