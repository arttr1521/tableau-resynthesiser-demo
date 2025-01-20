import random


def generate_random_matrix(n, m):
    """Generates a random boolean matrix with dimensions 2n x m."""
    return [[random.choice([True, False]) for _ in range(m)] for _ in range(2 * n)]

def print_changed_rows(matrix, changed_rows=[]):
    """Pretty prints the matrix with optional changed rows highlighted."""
    for i, row in enumerate(matrix):
        if i == len(matrix) // 2:
            print('-' * (len(row) * 2 - 1))  # Barrier line in the middle
        if i in changed_rows:
            print(' '.join(['1' if val else '0' for val in row]) + '  <--')
        else:
            print(' '.join(['1' if val else '0' for val in row]))
    print()

def main():
    # Fixing the random seed for reproducibility
    random.seed(42)

    # Step 1: Get N and M and generate the matrix
    try:
        n = int(input("Enter the number of rows (N): "))
        m = int(input("Enter the number of columns (M): "))
    except ValueError:
        print("Please enter valid integers for N and M.")
        return

    matrix = generate_random_matrix(n, m)

    print("Initial Matrix:")
    print_changed_rows(matrix)

    # Step 2 & 3: Accept operations
    operation_log = []
    while True:
        command = input("Enter operation or 'exit': ").strip().lower()
        changed_rows = []
        if command == 'exit':
            print("Exiting program.")
            print("Operations log:")
            for operation in operation_log:
                print(operation)
            break
        elif command == 'history':
            print("Operations log:")
            for operation in operation_log:
                print(operation)
        else:
            operation_log.append(command)
            try:
                parts = command.split()
                if parts[0] == 'swap' and len(parts) == 3:
                    row1, row2 = int(parts[1]), int(parts[2])
                    if 0 <= row1 < 2 * n and 0 <= row2 < 2 * n:
                        matrix[row1], matrix[row2] = matrix[row2], matrix[row1]
                        print(f"Rows {row1} and {row2} swapped.")
                        changed_rows = [row1, row2]  # Highlight both swapped rows
                    else:
                        print("Invalid row indices.")
                elif parts[0] == 'reverse' and len(parts) == 2:
                    row = int(parts[1])
                    if 0 <= row < 2 * n:
                        matrix[row] = matrix[row][::-1]
                        print(f"Row {row} reversed.")
                        changed_rows = [row]
                    else:
                        print("Invalid row index.")
                elif parts[0] == 'set' and len(parts) == 3:
                    row = int(parts[1])
                    value = bool(int(parts[2]))
                    if 0 <= row < 2 * n:
                        matrix[row] = [value] * m
                        print(f"Row {row} set to {'1' if value else '0' }.")
                        changed_rows = [row]
                    else:
                        print("Invalid row index.")
                elif parts[0] == 'add' and len(parts) == 3:
                    row1, row2 = int(parts[1]), int(parts[2])
                    if 0 <= row1 < 2 * n and 0 <= row2 < 2 * n:
                        matrix[row2] = [(val1 ^ val2) for val1, val2 in zip(matrix[row1], matrix[row2])]
                        print(f"Row {row1} XORed with row {row2}.")
                        changed_rows = [row1, row2]
                    else:
                        print("Invalid row indices.")
                elif parts[0] == 's' and len(parts) == 2:
                    row = int(parts[1])
                    if 0 <= row < n:
                        matrix[row] = [(val1 ^ val2) for val1, val2 in zip(matrix[row], matrix[row + n])]
                        print(f"Row {row} added with row {row + n}.")
                        changed_rows = [row]
                    else:
                        print("Invalid row index.")
                elif parts[0] == 'h' and len(parts) == 2:
                    row = int(parts[1])
                    if 0 <= row < n:
                        matrix[row], matrix[row + n] = matrix[row + n], matrix[row]
                        print(f"Rows {row} and {row + n} swapped.")
                        changed_rows = [row, row + n]
                    else:
                        print("Invalid row index.")
                elif parts[0] == 'cx' and len(parts) == 3:
                    row1, row2 = int(parts[1]), int(parts[2])
                    if 0 <= row1 < n and 0 <= row2 < n:
                        matrix[row1] = [(val1 ^ val2) for val1, val2 in zip(matrix[row1], matrix[row2])]
                        matrix[row2 + n] = [(val1 ^ val2) for val1, val2 in zip(matrix[row1 + n], matrix[row2 + n])]
                        print(f"Controlled-X operation between rows {row1} and {row2} performed.")
                        changed_rows = [row1, row2 + n]
                    else:
                        print("Invalid row indices.")
                else:
                    print("Invalid command. Supported operations: 'swap', 'reverse', 'set', 'add', 's', 'h', 'cx'.")
            except (ValueError, IndexError):
                print("Invalid command format.")
        
        # Print the updated matrix with the changed rows highlighted
        print_changed_rows(matrix, changed_rows)

if __name__ == "__main__":
    main()
