from unittest import TestCase
from tictactoe import *


class Test(TestCase):

    def test_generates_successor(self):
        before = 'XOX......'
        after = successor(before, 'O', 5)
        self.assertEqual('XOX..O...', after)

    def test_generates_legal_moves(self):
        before = '.OXO...X.'
        moves = legal_moves(before, 'X')
        self.assertEqual((0, 4, 5, 6, 8), moves)

    def test_finds_win_for_X(self):
        board = 'XXXO.O.O.'
        result = winner(board)
        self.assertEqual(1, result)

    def test_finds_value_of_end_of_game(self):
        board = 'XO.OX...X'
        result = value(board, 'O')
        self.assertEqual(1, result)

    def test_looks_ahead_one_move(self):
        board = 'XX.OO..X.'
        result = value(board, 'O')
        self.assertEqual(-1, result)

    def test_looks_ahead_multiple_moves(self):
        board = '.XX.O.OX.'
        result = value(board, 'O')
        self.assertEqual(-1, result)

    def test_does_not_find_legal_moves_when_game_is_over(self):
        board = 'XXXO.O.O.'
        result = legal_moves(board, 'X')
        self.assertEqual((), result)

    def test_finds_best_move(self):
        board = 'XOX.OO.X.'
        result = best_move(board, 'X')
        self.assertEqual(3, result)