import pygame
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode
import random
from datetime import datetime
from universe import Universe
from heavenly_body import *

# set up the colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# constants
G = 6.674e-11 # N kg-2 m^2
# Earth_Mass = 5.972e24 # kg
# Moon_Mass = 7.34767309e22 # kg
Distance = 384400000. # m
Max_Speed = 300

def main():
    print('Press q to quit')

    random.seed(0)
    
    # Initializing pygame
    pygame.init()
    win_width = 1500
    win_height = 800
    dt = 100
    paused = False

    screen = pygame.display.set_mode((win_width, win_height))  # Top left corner is (0,0)
    pygame.display.set_caption('Heavenly Bodies')

    mass, radius, temp, luminosity = HeavenlyBody.generate_star()
    star = Star([win_width/2/scale, win_height/2/scale], [0,0], mass, radius, temp, luminosity)
    universe = Universe(star, [])
    num_planets = 9
    # star_distance = star.pos[0] + star.radius
    star_distance = star.radius
    print(star.pos)

    for i in range(num_planets):
        star_distance += random.uniform(50000, 100000)

        if i >= num_planets *2/3:
            planet_type = "ice_giant"
        elif i >= num_planets / 3:
            planet_type = "gas_giant"
        else:
            planet_type = "rocky"


        if i % 4 == 0:
            pos = [star.pos[0] + star_distance, star.pos[1]]
            dir = 1
        elif i % 4 == 1:
            pos = [star.pos[0], star.pos[1] + star_distance]
            dir = -1
        elif i % 4 == 2:
            pos = [star.pos[0] - star_distance, star.pos[1]]
            dir = 1
        else:
            pos = [star.pos[0], star.pos[1] - star_distance]
            dir = -1

        #pos = [star.pos[0] + star_distance, star.pos[1]]
        vel = np.array([1e4*(-1**i), 0])
        
        mass, radius, temp, albedo = HeavenlyBody.generate_planet(planet_type)
        
        if i % 2 == 0:
            pos[0] += dir*radius
        else:
            pos[1] += dir*radius

        universe.planets.append(Planet(pos, vel, mass, radius, temp, albedo, planet_type))
        universe.planets[i].vel = Universe.calculate_orbital_velocity(star, universe.planets[i], star_distance)
        # universe.planets[i].vel = np.array([1.0e4 * i+1, 0.0])
        # universe.planets[i].vel = vel
        print(universe.planets[i].pos)
        print(universe.planets[i].vel)
        star_distance += 2*radius

    total_frames = 1000000
    iter_per_frame = 1

    frame = 0
    while frame < total_frames:
        if False:
            print ('Frame number', frame)        

        event = pygame.event.poll()
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit(0)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_q:
            pygame.quit()
            sys.exit(0)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_p:
            paused = not paused
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_UP:
            universe.add_asteroid("top", win_width, win_height)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_LEFT:
            universe.add_asteroid("left", win_width, win_height)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_RIGHT:
            universe.add_asteroid("right", win_width, win_height)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_DOWN:
            universe.add_asteroid("bottom", win_width, win_height)
        else:
            pass

        if not paused:
            universe.update()
            if frame % 10 == 0:
                universe.update_trails()
            if frame % 200 == 0:
                for planet in universe.planets:
                    print(planet.type, planet.temp)
            if frame % iter_per_frame == 0:
                screen.fill((20,0,40)) # clear the background
                universe.draw(screen)
                pygame.display.flip()
            frame += 1

    pygame.quit()

if __name__ == '__main__':
    main(),