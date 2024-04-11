import argparse
import numpy as np


def main():
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
    parser.add_argument(
        '--coreg',
        '-c',
        type=str,
        help=
        'Coregistration folder path. If provided, will warp to MNI space and return AAL regions.'
    )

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
    print('CT space:')
    print('   x' + ' ' * 10 + 'y' + ' ' * 10 + 'z')
    print('-' * 10 + ' ' + '-' * 10 + ' ' + '-' * 10)
    for i, pt in enumerate(trajectory):
        print(f'{pt[0]: 10.5f} {pt[1]: 10.5f} {pt[2]: 10.5f}')

    # print trajectory in template space
    if args.coreg:
        print('\n' * 3 + '... warping to MNI space ...' + '\n' * 2)
        from .coreg import warpcoords_ct_to_MNI
        from .fuzzyquery_aal import lookup_aal_region

        # transform trajectory to template space
        mni_coords = warpcoords_ct_to_MNI(np.array(trajectory), args.coreg)

        print('MNI Template space:')
        print('   x' + ' ' * 10 + 'y' + ' ' * 10 + 'z' + ' ' * 7 +
              'AAL region')
        print('-' * 10 + ' ' + '-' * 10 + ' ' + '-' * 10 + ' ' + '-' * 10)
        for i, pt in enumerate(mni_coords):
            aal = lookup_aal_region(pt.tolist(), fuzzy_dist=10)[1]
            print(f'{pt[0]: 10.5f} {pt[1]: 10.5f} {pt[2]: 10.5f} {aal}')


if __name__ == '__main__':
    main()
