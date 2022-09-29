# @author Skye Russ

"""
If you've already placed all n queens, return the solution
Otherwise, for each row in the next column:
If you can solve the puzzle by placing the next queen in that row and then recurring, return that solution
If no solution was found, return False
"""


def check_space(r: int, c: int, board: tuple):
    """
    checks if the space passed is a valid queen location
    :param r: row to be filled
    :param c: column to be filled
    :param board: list of current queens rows (in column order 0 -> n)
    :return: False if not valid, True if valid
    """
    if r in board:  # if row already taken
        return False
    i = 1  # start at 1 less than the current row
    while i <= c:
        if c - i >= 0 and (board[c - i] == r - i or board[c - i] == r + i):  # check if column before available and if queen in diagonals
            return False
        i += 1
    return True


def find_solution(limit, c, queen_list):
    """
    Recursively finds if there is a solution at given column with a limit of columns
    :param limit: board size (column limit)
    :param c: current column number
    :param queen_list: rows of current queens placed
    :return: False if no possibilities, otherwise return tuple with list of rows for queens to fill
    """
    for r in range(limit):
        if check_space(r, c, queen_list):
            if limit == c+1:  # if on last row
                temp = queen_list + (r,)
            else:
                temp = find_solution(limit, c + 1, (queen_list + (r,)))
            if not temp:  # if recursive call comes back false, try the next row
                continue
            else:
                return temp
    return False


def nqueens(num):
    """
    returns tuple of rows for queens to be in with board size num (False if no options)
    :param num: number of queens / size of board
    :return: tuple of rows for queens or False if no possible answers
    """
    return find_solution(num, 0, ())
