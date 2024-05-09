import matplotlib.pyplot as plt
import numpy as np

def draw_upper_triangular_matrix_with_triangular_sylvester(n, square_size, step=0):
    fig, ax = plt.subplots(figsize=(5, 5))
    
    for j in range(n-step):
        for i in range(j, n-step):
            square = plt.Rectangle((i, j), square_size, square_size, fc='#CCCCCC', ec='black')
            ax.add_patch(square)

    for i in range(n-step):
        square = plt.Rectangle((n-step, i), square_size, square_size, fc='blue', ec='black')
        ax.add_patch(square)

    for i in range(n-step, n):
        square = plt.Rectangle((i, i), square_size, square_size, fc='yellow', ec='black')
        ax.add_patch(square)

    square = plt.Rectangle((n-step, n-step), square_size, square_size, fc='red', ec='black')
    ax.add_patch(square)

    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(f'assets/sylvester_{n}_{step}.pdf', bbox_inches='tight')
    plt.close()

def draw_sylvester_solve(n, square_size):
    fig, ax = plt.subplots(figsize=(5, 5))
    
    for j in range(n-1):
        for i in range(j, n-1):
            square = plt.Rectangle((i, j), square_size, square_size, fc='#CCCCCC', ec='black')
            ax.add_patch(square)
    
    ax.set_xlim(0, n-1)
    ax.set_ylim(0, n-1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(f'assets/sylvester_solve_{n}_lhs.pdf', bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 5))

    for i in range(n-1):
        square = plt.Rectangle((0, i), square_size, square_size, fc='blue', ec='black')
        ax.add_patch(square)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, n-1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'assets/sylvester_solve_{n}_rhs_u.pdf', bbox_inches='tight')
    plt.close()


    fig, ax = plt.subplots(figsize=(5, 5))

    for i in range(n-1):
        square = plt.Rectangle((0, i), square_size, square_size, fc='green', ec='black')
        ax.add_patch(square)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, n-1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'assets/sylvester_solve_{n}_solution.pdf', bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 5))

    square = plt.Rectangle((0, 0), square_size, square_size, fc='red', ec='black')
    ax.add_patch(square)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(f'assets/sylvester_solve_{n}_rhs_l.pdf', bbox_inches='tight')
    plt.close()
    
    

# Usage
n = 10
square_size = 1
# for step in range(0, n+1):
#     draw_upper_triangular_matrix_with_triangular_sylvester(n, square_size, step)

draw_sylvester_solve(n, square_size)