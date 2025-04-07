from heavenly_body import *
import numpy as np

class Collision:
    @staticmethod
    def range(universe, planet):
        pass

    @staticmethod
    def check_all(universe, body):
        colliding = []
        
        if Collision.check(universe.star, body):
            colliding.append(universe.star)
        
        for planet in universe.planets:
            if planet == body:
                continue

            elif Collision.check(planet, body):
                colliding.append(planet)
        return colliding

    @staticmethod
    def check(body1, body2):
        distance = np.linalg.norm(body1.pos - body2.pos)
        return distance <= body1.radius + body2.radius

    @staticmethod
    def response(body1, body2):
        # Vector between the two bodies
        delta_pos = body2.pos - body1.pos
        distance = np.linalg.norm(delta_pos)

        # If the bodies are colliding, we need to calculate the response
        if distance <= body1.radius + body2.radius:
            # Normalize the vector to get the direction of the collision
            collision_normal = delta_pos / distance

            # Calculate relative velocity along the collision normal
            relative_velocity = body2.vel - body1.vel
            velocity_along_normal = np.dot(relative_velocity, collision_normal)

            # If the bodies are moving towards each other, we need to calculate new velocities
            if velocity_along_normal < 0:  # Only handle if they are moving towards each other
                # Coefficients of restitution (elastic collision: e = 1, inelastic < 1)
                restitution = 1.0  # Assuming perfectly elastic collision

                # Calculate impulse scalar
                impulse_magnitude = -(1 + restitution) * velocity_along_normal
                impulse_magnitude /= (1 / body1.mass + 1 / body2.mass)

                overlap = body1.radius + body2.radius - distance
                impulse = impulse_magnitude * collision_normal
                return impulse, (overlap / 2 * collision_normal)