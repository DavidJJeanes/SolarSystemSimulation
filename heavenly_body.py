import random
import pygame
import numpy as np
from datetime import datetime

scale = 1/2000
SIGMA = 5.67e-8  # Stefan-Boltzmann constant (W/m^2 K^4)
# SIGMA = 5.0e-11
MAX_TEMP = 1e6  # 1 million K (as an upper limit)

class HeavenlyBody:
    def __init__(self, pos, vel, mass, radius):
        self.mass = mass
        self.radius = radius
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.grid_pos = 0

        #self.image = pygame.Surface([radius, radius])
        #self.image.fill((0,0,0))
        #pygame.draw.circle(self.image, (1, 1, 1), (radius, radius), radius, radius)

    def draw(self, screen):
        """Draw the body at its current position on the screen."""
        # Draw the image at the current position, adjusting for scaling
        pygame.draw.circle(screen, self.colour, (int(self.pos[0] * scale), int(self.pos[1] * scale)), int(self.radius*scale))

    @staticmethod
    def generate_planet(planet_type):
        if planet_type == "rocky":
            # Scaled similar to Earth
            radius = random.uniform(7000, 9000)
            mass = 5.97e24 * (radius / 6371) ** 3
            temp = random.uniform(250, 310)
            albedo = random.uniform(0.2, 0.3)
        elif planet_type == "gas_giant":
            # Scaled similar to Jupiter
            radius = random.uniform(50000, 70000)
            mass = 1.9e27 * (radius / 69911) ** 3
            temp = random.uniform(90, 150)
            albedo = random.uniform(0.5, 0.6)
        elif planet_type == "ice_giant":
            # Scaled similar to Neptune
            radius = random.uniform(20000, 30000)
            mass = 1.0e26 * (radius / 24622) ** 3
            temp = random.uniform(50, 100)
            albedo = random.uniform(0.6, 0.7)
        return mass, radius, temp, albedo

    @staticmethod
    def generate_star():
        #mass = random.uniform(1.5e24, 2.5e24)
        mass = random.uniform(5.5e23, 7.5e23)
        radius = random.uniform(150000, 200000)
        temp = random.uniform(5000, 7000)
        # luminosity = random.uniform(0.8, 1.2)
        # luminosity = random.uniform(3.5e26, 4.5e26)
        luminosity = 1e26
        luminosity = 4 * np.pi * (radius * 1000) ** 2 * SIGMA * temp ** 4
        return mass, radius, temp, luminosity

    @staticmethod
    def generate_asteroid():
        pass


class Star(HeavenlyBody):
    def __init__(self, pos, vel, mass, radius, temp, luminosity):
        super().__init__(pos, vel, mass, radius)
        self.temp = temp
        self.luminosity = luminosity
        self.colour = (random.randrange(155, 255), random.randrange(155, 255), random.randrange(100))

class Planet(HeavenlyBody):
    def __init__(self, pos, vel, mass, radius, temp, albedo, type):
        super().__init__(pos, vel, mass, radius)
        self.temp = temp
        self.albedo = albedo
        self.type = type
        
        if type == "rocky":
            self.colour = (random.randrange(100, 255), random.randrange(30, 155), random.randrange(30, 155))
            self.heat_capacity = random.uniform(750, 850)
        elif type == "gas_giant":
            self.colour = (random.randrange(30, 155), random.randrange(100, 255), random.randrange(30, 155))
            self.heat_capacity = random.uniform(1100, 1300)
        else:
            self.colour = (random.randrange(30, 155), random.randrange(30, 155), random.randrange(100, 255))
            self.heat_capacity = random.uniform(1800, 2200)

        self.trail = [self.pos] * 100

    def draw_trail(self, screen):
        fade = 0.2
        for trail in self.trail:
                fade += 0.005
                faded_colour = tuple(min(int(c * fade), 255) for c in self.colour) 
                pygame.draw.circle(screen, faded_colour, (int(trail[0] * scale), int(trail[1] * scale)), int(self.radius*scale))

    """
    def update_temp(self, star, dt):
        distance = np.linalg.norm(self.pos - star.pos) * 1000
        # Calculate the incoming energy from the star using the inverse square law
        # Luminosity of the star is L_star (you can set a value for L_star)
        L_star = star.luminosity  # Luminosity of the star (W)
        flux = L_star / (4 * np.pi * distance**2)  # Flux received by the planet
        print(f"Star Luminosity: {star.luminosity} W")
        print(f"Distance to star: {distance} meters")
        print(f"Planet Radius: {self.radius} km")
        print(f"Planet Albedo: {self.albedo}")
        print(f"Flux: {flux} W/m^2")

        # Calculate energy absorbed (absorbed by the planet's surface)
        absorbed_energy = (1 - self.albedo) * flux * np.pi * ((self.radius * 1000)**2) / 1e7 # Total energy absorbed
        print(f"Absorbed Energy: {absorbed_energy} J")

        # Calculate energy radiated by the planet (blackbody radiation)
        radiated_energy = SIGMA * (self.radius * 1000)**2 * 4 * np.pi * (self.temp ** 4)  # Stefan-Boltzmann law

        # Rate of change of energy (stored energy) due to incoming energy and radiation
        delta_energy = absorbed_energy - radiated_energy  # Net energy gain or loss
        print(f"Radiated Energy: {radiated_energy} J")

        # Update the planet's energy based on the net energy and thermal inertia (heat capacity)
        self.temp += delta_energy * dt / self.heat_capacity  # Update stored energy
        self.temp = min(self.temp, MAX_TEMP)  # Clamp temperature
    """
    
    def update_temp(self, star, dt):
        # Equilibrium temp
        distance = np.linalg.norm(self.pos - star.pos)**1.04 * 1000 * 80
        
        self.temp = (((1 - self.albedo) * star.luminosity) / (16*np.pi*SIGMA*(distance**2)))**0.25
        # print(self.type, self.temp)

class Asteroid(HeavenlyBody):
    def __init__(self, pos, vel, mass, radius):
        super().__init__(pos, vel, mass, radius)
        # self.colour = (random.randrange(255), random.randrange(255), random.randrange(255))
        self.colour = (255, 255, 255)
        self.trail = [self.pos] * 50

    def draw_trail(self, screen):
        fade = 0.1
        for trail in self.trail:
                fade += 0.02
                faded_colour = tuple(min(int(c * fade)+20, 255) for c in self.colour) 
                pygame.draw.circle(screen, faded_colour, (int(trail[0] * scale), int(trail[1] * scale)), int(self.radius*scale))