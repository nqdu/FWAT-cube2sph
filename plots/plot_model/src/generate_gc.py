import numpy as np
import sys

def generate_gc(lon1, lat1, lon2, lat2, num_points=100, radius=6371.0):
    """
    Generate a great-circle path between two points using only NumPy.

    Parameters:
    - lon1, lat1: coordinates of the first point (in degrees)
    - lon2, lat2: coordinates of the second point (in degrees)
    - num_points: number of interpolated points
    - radius: sphere radius in meters (default: Earth's mean radius)

    Returns:
    - lons: longitudes along the path (degrees)
    - lats: latitudes along the path (degrees)
    - dists: distances from the first point (meters)
    """
    # Convert to radians
    lon1, lat1 = np.radians([lon1, lat1])
    lon2, lat2 = np.radians([lon2, lat2])

    # Convert to Cartesian coordinates
    def sph2cart(lon, lat):
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat)
        ])

    A = sph2cart(lon1, lat1)
    B = sph2cart(lon2, lat2)

    # Angle between A and B
    omega = np.arccos(np.clip(np.dot(A, B), -1.0, 1.0))

    if omega == 0:
        # Points are the same
        lats = np.full(num_points, np.degrees(lat1))
        lons = np.full(num_points, np.degrees(lon1))
        dists = np.zeros(num_points)
        return lons, lats, dists

    # Interpolation fractions
    f = np.linspace(0, 1, num_points)

    # Slerp (spherical linear interpolation)
    sin_omega = np.sin(omega)
    points = (np.sin((1 - f) * omega)[:, None] * A + np.sin(f * omega)[:, None] * B) / sin_omega

    # Convert back to lat/lon
    lats = np.degrees(np.arcsin(points[:, 2]))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))

    # Distance from first point along the arc
    dists = f * omega * radius

    return lons, lats, dists


def main():
    if len(sys.argv)!=7:
        print("Usage ./this lonmin lonmax latmin latmax n outfile")
        exit(1)
    
    lonmin,lonmax,latmin,latmax = map(lambda x:float(x),sys.argv[1:5])
    n = int(sys.argv[5])
    outfile = sys.argv[6]
    data = np.zeros((n,3))

    data[:,0],data[:,1],data[:,2] = generate_gc(lonmin,latmin,lonmax,latmax,n)

    np.savetxt(outfile,data,fmt='%f')

if __name__ == "__main__":
    main()
