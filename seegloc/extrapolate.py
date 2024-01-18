import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description=
    'Trajectory extrapolation tool: Use to manually recover trajectories from CT images when individual electrodes are not descernable.'
)
parser.add_argument('spacing',
                    type=float,
                    help='Spacing between electrodes in mm')
parser.add_argument('n_electrodes', type=int, help='Number of electrodes')
parser.add_argument('start',
                    type=float,
                    nargs=3,
                    help='Electrode tip in mm (CT space)')
parser.add_argument('point',
                    type=float,
                    nargs=3,
                    help='Any other point on the trajectory (CT space)')

args = parser.parse_args()

# get norm direction from start to point
start = np.array(args.start)
point = np.array(args.point)

direction = point - start
direction /= np.linalg.norm(direction)

# extrapolate trajectory from start
trajectory = [start]
for i in range(args.n_electrodes - 1):
    trajectory.append(trajectory[-1] + direction * args.spacing)

# print trajectory
print('   x' + ' ' * 10 + 'y' + ' ' * 10 + 'z')
print('-' * 10 + ' ' + '-' * 10 + ' ' + '-' * 10)
for i, pt in enumerate(trajectory):
    print(f'{pt[0]: 10.5f} {pt[1]: 10.5f} {pt[2]: 10.5f}')
