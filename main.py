from synthesiser import tableau_resynthesis, Btor2generator
import argparse
import random

def read_tableau(filename):
    """
    Reads a tableau from a file and returns it as a list of lists.
    The input file format is expected to be:
      <number of rows> <number of columns>
      row1
      row2
      ...

    Parameters:
        filename (str): The path to the tableau file.

    Returns:
        List[List[int]]: The tableau as a list of lists of integers.
    """
    with open(filename, 'r') as file:
        # Read the first line to get the number of rows and columns
        header = file.readline().strip()
        if not header:
            raise ValueError('The file is empty or header is missing.')
        try:
            num_rows, num_cols = map(int, header.split())
        except ValueError:
            raise ValueError('Header must contain two integers: <number of rows> <number of columns>')

        tableau = []
        for i in range(num_rows):
            line = file.readline().strip()
            if len(line) < num_cols:
                raise ValueError(f'Line {i+2} does not contain enough characters for {num_cols} columns.')
            # Convert each character in the line to int
            row = [int(c) for c in line[:num_cols]]
            tableau.append(row)
    return tableau

def write_tableau(tableau, filename):
    """
    Writes the tableau to a file.
    The output file format will be:
      <number of rows> <number of columns>
      row1 (as a string of digits)
      row2
      ...

    Parameters:
        tableau (List[List[int]]): The tableau to write.
        filename (str): The output file name.
    """
    num_rows = len(tableau)
    num_cols = len(tableau[0]) if tableau else 0
    with open(filename, 'w') as file:
        file.write(f"{num_rows} {num_cols}\n")
        for row in tableau:
            file.write("".join(map(str, row)) + "\n")

def generate_random_tableau(num_qubits, num_rotations, seed=None):
    """
    Generates a random tableau for testing with an optional fixed random seed,
    ensuring that no column is entirely zero.
    Generates a random tableau for testing with an optional fixed random seed,
    ensuring that no column is entirely zero.

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

    # Ensure that no column is all zeros
    for col in range(cols):
        if all(tableau[row][col] == 0 for row in range(rows)):  
            # Pick a random row and set it to 1
            random_row = random.randint(0, rows - 1)
            tableau[random_row][col] = 1

    # Ensure that no column is all zeros
    for col in range(cols):
        if all(tableau[row][col] == 0 for row in range(rows)):  
            # Pick a random row and set it to 1
            random_row = random.randint(0, rows - 1)
            tableau[random_row][col] = 1

    return tableau

def run_ubmc(tableau, max_sat):
    print("Running UBMC mode")
    # Initialize the tableau_resynthesis object
    resynthesizer = tableau_resynthesis(tableau)
    # resynthesizer.setMaxDepth(10)
    resynthesizer.UBMC()
    # If UBMC succeeds, print success message
    print("\nUBMC completed successfully.")

    resynthesizer.print_graph()
    
    # Print variables for debugging
    print("\nGenerated Variables:")
    resynthesizer.print_variables()

    # Print the CNF clauses for debugging
    print("\nGenerated CNF Clauses:")
    resynthesizer.print_clauses(detail=True)
    
    if max_sat:
        resynthesizer.optimize_result()
        resynthesizer.print_detail(style="maxsat")
    else:
        # resynthesizer.print_result(style="cex")
        resynthesizer.print_result(style="detail")


def run_bmc(tableau, depth, max_sat):
    print(f"Running BMC mode with depth {depth}")
    
    # Initialize the tableau_resynthesis object
    resynthesizer = tableau_resynthesis(tableau)
    resynthesizer.BMC(depth)
    # If BMC succeeds, print success message
    print("\nBMC completed successfully.")
    # resynthesizer.print_result(style="cex")
    resynthesizer.print_result(style="detail")

    resynthesizer.print_graph()

    # Print variables for debugging
    print("\nGenerated Variables:")
    resynthesizer.print_variables()

    # Print the CNF clauses for debugging
    print("\nGenerated CNF Clauses:")
    resynthesizer.print_clauses(detail=True)

    if max_sat:
        resynthesizer.optimize_result()
        resynthesizer.print_detail(style="maxsat")
    else:
        resynthesizer.print_result(style="detail")



def run_cnf(tableau, depth, output):
    print(f"Running CNF mode with output file {output}")
    # To be implemented
    resynthesizer = tableau_resynthesis(tableau)
    resynthesizer.build_bmc_cnf(depth, output)

def run_btor2(tableau, output):
    print(f"Running BTOR2 mode with output file {output}")
    resynthesizer = tableau_resynthesis(tableau)
    resynthesizer.generate_btor2_file(output)
    resynthesizer.print_graph()

def main():
    parser = argparse.ArgumentParser(description="Tableau Resynthesiser Demo")
    parser.add_argument('mode', choices=['ubmc', 'bmc', 'cnf', 'btor2'], help='Select mode: ubmc, bmc, cnf, btor2')
    parser.add_argument('--depth', type=int, default=10, help='Depth parameter required for BMC mode (only applicable to bmc mode)')
    parser.add_argument('--output', type=str, help='Output file for CNF or BTOR2 modes (only applicable to cnf and btor2 modes)')
    parser.add_argument('--seed', type=int, default=None, help='Seed for random number generator')
    parser.add_argument('--num_qubits', type=int, default=4, help='Number of qubits for random generation')
    parser.add_argument('--num_rotations', type=int, default=6, help='Number of rotations for random generation')
    parser.add_argument('--input', type=str, help='Input file for tableau, random for random generation')
    parser.add_argument('--max_sat', type=bool, default=False, help='Use max SAT solver')
    args = parser.parse_args()

    # Define a random tableau
    if args.input == "random":
        tableau = generate_random_tableau(args.num_qubits, args.num_rotations, seed=args.seed)
    else:
        tableau = read_tableau(args.input)

    print("Input tableau is: ")
    for row in tableau:
        print(row)

    if args.mode == 'ubmc':
        run_ubmc(tableau, args.max_sat)
    elif args.mode == 'bmc':
        run_bmc(tableau, args.depth, args.max_sat)
    elif args.mode == 'cnf':
        if not args.output:
            run_cnf(tableau, args.depth, "benchmark/cnf/output.cnf")
        else:
            run_cnf(tableau, args.depth, args.output)
    elif args.mode == 'btor2':
        if not args.output:
            run_btor2(tableau, "benchmark/btor2/output.btor2")
        else:
            run_btor2(tableau, args.output)

if __name__ == "__main__":
    main()