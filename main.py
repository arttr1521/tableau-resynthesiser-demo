from synthesiser import tableau_resynthesis

import random

def generate_random_tableau(num_qubits, num_rotations, seed=None):
    """
    Generates a random tableau for testing with an optional fixed random seed.

    Parameters:
        num_qubits (int): The number of qubits (N).
        num_rotations (int): The number of rotations (M).
        seed (int, optional): A seed for the random number generator. Default is None.

    Returns:
        List[List[int]]: A 2N x M random tableau with 0s and 1s.
    """
    if seed is not None:
        random.seed(seed)  # Set the random seed for reproducibility
    rows = 2 * num_qubits  # Tableau has 2N rows
    cols = num_rotations   # Tableau has M columns

    # Generate a 2N x M tableau with random 0s and 1s
    tableau = [[random.randint(0, 1) for _ in range(cols)] for _ in range(rows)]

    return tableau


def test_synthesiser(tableau):
    """
    Test function for the tableau_resynthesis class.
    Validates the UBMC process, constraints, and SAT solver setup.
    """

    # Initialize the tableau_resynthesis object
    resynthesizer = tableau_resynthesis(tableau)
    resynthesizer.setMaxDepth(10)
    resynthesizer.UBMC()

    resynthesizer.print_graph()

    # Print variables for debugging
    print("\nGenerated Variables:")
    resynthesizer.print_variables()

    # Print the CNF clauses for debugging
    print("\nGenerated CNF Clauses:")
    resynthesizer.print_clauses(detail=True)


    # If UBMC succeeds, print success message
    print("\nUBMC completed successfully.")
    resynthesizer.print_result(style="cex")
    resynthesizer.print_result(style="detail")

# Call the test function
if __name__ == "__main__":
    
    # # Define a simple tableau for testing
    # tableau = [
    #     [1, 0, 0, 1],
    #     [1, 0, 1, 0],
    #     [1, 1, 0, 0],
    #     [0, 1, 0, 1]
    # ]
    
    # Define a random tableau for testing
    tableau = generate_random_tableau(4, 4, seed=25)
    # print("Input tableau is: ")
    # for row in tableau:
    #     print(row)

    test_synthesiser(tableau)
