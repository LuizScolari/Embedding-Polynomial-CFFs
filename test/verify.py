def verify_uniform_sum_matrix(file_path, expected_row_sum, expected_col_sum):
    """
    Reads a binary matrix from a file and verifies if the sum of 1s in each row and
    column matches a single expected value for all rows and columns.

    Args:
        file_path (str): Path to the .txt file containing the matrix.
        expected_row_sum (int): Expected sum value for ALL rows.
        expected_col_sum (int): Expected sum value for ALL columns.
    """
    matrix = []
    # Read file and build matrix
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Convert each element in the line to int and add to matrix
                matrix.append([int(x) for x in line.strip().split()])
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return
    except ValueError:
        print("Error: File contains non-integer characters.")
        return

    if not matrix:
        print("Matrix is empty.")
        return

    num_rows = len(matrix)
    # Assume all rows have the same number of columns
    num_cols = len(matrix[0])

    print("--- Checking Rows ---")
    check = True
    # Iterate over each row to check the sum
    for i, row in enumerate(matrix):
        current_row_sum = sum(row)
        if current_row_sum != expected_row_sum:
            print(f"Row {i}: Expected: {expected_row_sum}, Found: {current_row_sum}")
            check = False
    if check == True:
        print("Valid!")
            
    check = True
    print("\n--- Checking Columns ---")
    # Iterate over each column to check the sum
    for j in range(num_cols):
        current_col_sum = sum(matrix[i][j] for i in range(num_rows))
        if current_col_sum != expected_col_sum:
            print(f"Column {j}: Expected: {expected_col_sum}, Found: {current_col_sum}")
            check = False
    if check == True:
        print("Valid!")


file_path = 'saida.txt'
expected_row_sum = 16
expected_col_sum = 16

print(f"Verifying matrix with expected sum of '{expected_row_sum}' for rows and '{expected_col_sum}' for columns.\n")

verify_uniform_sum_matrix(
    file_path,
    expected_row_sum,
    expected_col_sum
)