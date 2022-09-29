def successor(board, player, index):
    """
    Returns board if player played there or False if invalid move
    :param board: string board
    :param player: "X" or "O"
    :param index: where to put letter
    :return: string with successor move
    """
    return board[:index] + player + board[index + 1:]


def legal_moves(board, player):
    """
    returns possible legal moves for player
    :param board: string board
    :param player: "X" or "O"
    :return: sequence of possible legal moves player can take
    """
    if winner(board) != 0:
        return ()
    return tuple([i for i in range(len(board)) if board[i] == "."])


def winner(board):
    """
    Returns if the board has a winner, and who it is
    :param board: string current board state
    :return: 1 if "X" has won, -1 if "O" has won, or 0 otherwise
    """
    winning_lines = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6))
    for a, b, c in winning_lines:
        if board[a] == board[b] == board[c]:
            if board[a] == "X":
                return 1
            elif board[a] == "O":
                return -1
    return 0


def value(board, player):
    """
    Computes the value of the board for the parameter player
    :param board: string board
    :param player: "X" or "O"
    :return: -1 for O, 1 for X, 0 for tie
    """
    moves = legal_moves(board, player)
    if moves:
        if player == 'X':
            best = max
        else:
            best = min
        return best([value(successor(board, player, move), opposite(player)) for move in moves])
    return winner(board)


def best_move(board, player):
    """
    Computes the best move for the parameter player
    :param board: string board
    :param player: "X" or "O"
    :return: integer spot to fill
    """
    moves = legal_moves(board, player)
    val = 0
    best = -1
    for move in moves:
        temp = value(successor(board, player, move), player)
        if player == "X" and temp > val:
            best = move
            val = temp
        elif player == "O" and temp < val:
            best = move
            val = temp
    return best


def opposite(player):
    if player == "X":
        return "O"
    return "X"


def main(board, player):
    while legal_moves(board, player):
        print(f'{board[0:3]}\n{board[3:6]}\n{board[6:9]}\n')
        if player == "X":
            move = best_move(board, player)
        else:
            move = int(input("Your move: "))
        board = successor(board, player, move)
        player = opposite(player)
    print(f'{board[0:3]}\n{board[3:6]}\n{board[6:9]}\n')


if __name__ == '__main__':
    main(".........", "X")
