import numpy as np
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from ProjectAstro2 import planet_state_jd, JD_from_UTC

def get_phase_angle(planet_name, jd):
    r_earth, _ = planet_state_jd('Earth', jd, in_km=False)
    r_planet, _ = planet_state_jd(planet_name, jd, in_km=False)
    
    # Phase angle: angle between Sun-Earth and Sun-Planet vectors
    cos_alpha = np.dot(r_earth, r_planet) / (np.linalg.norm(r_earth) * np.linalg.norm(r_planet))
    angle = np.degrees(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))
    
    # Check if planet is "ahead" or "behind" Earth in longitude
    # (Simple 2D approximation in ecliptic plane)
    lon_earth = np.arctan2(r_earth[1], r_earth[0])
    lon_planet = np.arctan2(r_planet[1], r_planet[0])
    diff = (np.degrees(lon_planet - lon_earth)) % 360.0
    
    return angle, diff

def main():
    launch_date = (2026, 1, 5)
    jd = JD_from_UTC(*launch_date)
    planets = ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
    
    print(f"Planetary Alignment on {launch_date[2]:02d}/{launch_date[1]:02d}/{launch_date[0]}")
    print(f"{'Planet':<10} | {'Phase Angle (deg)':<20} | {'Relative Longitude (deg)':<25}")
    print("-" * 65)
    
    for planet in planets:
        angle, diff = get_phase_angle(planet, jd)
        print(f"{planet:<10} | {angle:<20.2f} | {diff:<25.2f}")

if __name__ == "__main__":
    main()
