from scipy.integrate import ode
import numpy as np
from collisions import *
import math

G = 6.674e-11 # N kg-2 m^2
Distance = 384400000.

win_width = 1500
win_height = 800

class Universe:
    def __init__(self, star, planets):
        self.w, self.h = 2.6*Distance, 2.6*Distance 
        self.star = star
        self.planets = planets
        self.dt = 0.02
        #self.dt = 0.01e6
        self.cur_time = 0.0
        self.asteroids = []
        

    def draw(self, screen):
        for asteroid in self.asteroids:
            asteroid.draw_trail(screen)
            asteroid.draw(screen)
        
        self.star.draw(screen)
        
        for planet in self.planets:
            planet.draw_trail(screen)
            planet.draw(screen)

    def update(self):
        self.cur_time += self.dt

        for planet in self.planets:
                d = [self.star.pos[0] - planet.pos[0], self.star.pos[1] - planet.pos[1]]
                r = np.linalg.norm(d)
                u = d / r
                f = u * G * planet.mass * self.star.mass / (r*r)

                accel = np.array(f) / planet.mass
                planet.pos = planet.pos + planet.vel * self.dt + 0.5 * accel * self.dt**2
                planet.vel = planet.vel + accel * self.dt  # Update the velocity

                planet.update_temp(self.star, self.dt)

                colliding = Collision.check_all(self, planet) 
                for other in colliding:
                    impulse = Collision.response(planet, other)
                
                    # If collision response returns a valid impulse, apply it
                    if impulse is not None:    
                        planet.vel -= impulse[0] / planet.mass
                        other.vel += impulse[0] / other.mass if not isinstance(other, Star) else 0
                        
                        planet.pos -= impulse[1]
                        other.pos += impulse[1] if not isinstance(other, Star) else 0

        for asteroid in self.asteroids:
                accel = np.array(f) / asteroid.mass
                asteroid.pos = asteroid.pos + asteroid.vel * self.dt + 0.5 * accel * self.dt**2

                colliding = Collision.check_all(self, asteroid) 
                for other in colliding:
                    if isinstance(other, Star) and (np.linalg.norm(asteroid.pos - other.pos) <= other.radius - (asteroid.radius*2))\
                                    or (math.fabs(asteroid.pos[0]) > (win_width / scale) + asteroid.radius)\
                                    or (math.fabs(asteroid.pos[1]) > (win_height / scale) + asteroid.radius):
                        self.asteroids.remove(asteroid)
                    elif not isinstance(other, Star):
                        impulse = Collision.response(asteroid, other)
                
                        # If collision response returns a valid impulse, apply it
                        if impulse is not None:    
                            asteroid.vel -= impulse[0] / asteroid.mass
                            other.vel += impulse[0] / other.mass
                        
                            asteroid.pos -= impulse[1]
                            other.pos += impulse[1]

    def update_trails(self):
        for asteroid in self.asteroids:
                for i in range(len(asteroid.trail) -1):
                    asteroid.trail[i] = asteroid.trail[i+1]
                asteroid.trail[-1] = asteroid.pos

        for planet in self.planets:
            for i in range(len(planet.trail) -1):
                planet.trail[i] = planet.trail[i+1]
            planet.trail[-1] = planet.pos


    def add_asteroid(self, direction, win_width, win_height):
        speed = 8e4
        if direction == "top":
            pos = [win_width/2 / scale, 0]
            vel = [0, speed]
        elif direction == "left":
            pos = [0, win_height/2 / scale]
            vel = [speed, 0]
        elif direction == "right":
            pos = [win_width / scale, win_height/2 / scale]
            vel = [-speed, 0]
        else:
            pos = [win_width/2 / scale, win_height / scale]
            vel = [0, -speed]
        self.asteroids.append(Asteroid(np.array(pos), np.array(vel), 5.0e25, 25000))

    """
    @staticmethod
    def calculate_orbital_velocity(star, planet):
        # Calculate the distance from the star
        r = np.linalg.norm(planet.pos - star.pos)
    
        # Calculate the orbital velocity using v = sqrt(G * M / r)
        velocity_magnitude = np.sqrt(G * star.mass / r)
    
        # Calculate the velocity vector perpendicular to the radius vector
        # Assuming a circular orbit, the velocity vector will be tangential
        direction = np.array([- (planet.pos[1] - star.pos[1]), planet.pos[0] - star.pos[0]])  # Perpendicular to radial direction
        direction = direction / np.linalg.norm(direction)  # Normalize the direction
    
        # Set the velocity of the planet
        return velocity_magnitude * direction
    """

    @staticmethod
    def calculate_orbital_velocity(star, planet, star_distance, eccentricity=random.uniform(0.4, 0.8)):
        # Calculate the distance from the star (r)
        semi_major_axis = star_distance*1.25
        r = np.linalg.norm(planet.pos - star.pos)
        
        # Semi-major axis (a) and eccentricity (e) are used to calculate the orbital velocity
        a = semi_major_axis  # Semi-major axis of the ellipse
        e = eccentricity  # Eccentricity of the orbit

        # Calculate the orbital velocity using v = sqrt(G * M * (2/r - 1/a))
        velocity_magnitude = np.sqrt(G * star.mass * (2/r - 1/a))
        
        # Now calculate the direction of the velocity.
        # In elliptical orbits, the velocity vector is not just perpendicular but depends on the position.
        
        # Calculate the unit vector along the radial direction
        radial_vector = planet.pos - star.pos
        radial_unit_vector = radial_vector / np.linalg.norm(radial_vector)
        
        # The velocity is always perpendicular to the radius vector in ideal orbital mechanics
        # For simplicity, assume a 2D orbit in the x-y plane:
        # Use a perpendicular vector in 2D (assuming the plane is not tilted):
        direction = np.array([-radial_unit_vector[1], radial_unit_vector[0]])  # Perpendicular to radial direction
        direction = direction / np.linalg.norm(direction) 

        velocity = velocity_magnitude * direction
        
        # Return the velocity as a vector
        return velocity

    """
    def calculate_orbital_velocity(star, planet, star_distance, eccentricity=random.uniform(0.4, 0.8), theta=np.pi/3):
        Calculate the orbital velocity of a planet in an elliptical orbit using the Vis-Viva equation.
        
        Parameters:
        - star: The star object containing the mass of the star (star.mass).
        - planet: The planet object containing its initial position (planet.pos) and velocity (planet.vel).
        - star_distance: The distance from the star to the planet at a reference position (in meters).
        - eccentricity: The eccentricity of the orbit (0 <= e < 1).
        - theta: The true anomaly (angle) at which the planet is located in radians.
        
        Updates the planet's velocity.

        # Gravitational constant
        G = 6.67430e-11
        
        # Calculate the semi-major axis (a)
        semi_major_axis = star_distance * 1.25
        a = semi_major_axis  # Semi-major axis of the ellipse
        e = eccentricity  # Eccentricity of the orbit
        
        # Calculate the radial distance (r) using the true anomaly and the eccentricity
        r = a * (1 - e**2) / (1 + e * np.cos(theta))
        
        # Calculate the orbital velocity using the Vis-Viva equation
        velocity_magnitude = np.sqrt(G * star.mass * (2 / r - 1 / a))
        
        # Now calculate the direction of the velocity.
        # First, calculate the radial vector from the planet to the star.
        radial_vector = planet.pos - star.pos
        radial_unit_vector = radial_vector / np.linalg.norm(radial_vector)
        
        # The velocity vector is perpendicular to the radial vector in the orbital plane.
        # In 2D (x-y plane), we can use a simple perpendicular direction:
        direction = np.array([-radial_unit_vector[1], radial_unit_vector[0]])  # Perpendicular to the radial direction
        direction = direction / np.linalg.norm(direction)  # Normalize the direction
    
        # Set the planet's velocity in the correct direction.
        return velocity_magnitude * direction
    """