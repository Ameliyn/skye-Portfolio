from unittest import TestCase
from othello import *


class Test(TestCase):

    def test_finds_flips(self):
        board = ('........',
                 '........',
                 '...OX...',
                 '..XXXX..',
                 '..XXXOXX',
                 '...O.XO.',
                 '........',
                 '........')
        result = flips(board, 'O', (4, 1))
        self.assertEqual({(3, 2), (4, 2), (4, 3), (4, 4)}, set(result))

    def test_finds_successor(self):
        board = ('........',
                 '........',
                 '...OX...',
                 '..XXXX..',
                 '..XXXOXX',
                 '...O.XO.',
                 '........',
                 '........')
        result = successor(board, 'O', (4, 1))
        correct = ('........',
                   '........',
                   '...OX...',
                   '..OXXX..',
                   '.OOOOOXX',
                   '...O.XO.',
                   '........',
                   '........')
        self.assertEqual(correct, result)

    def test_finds_legal_moves(self):
        board = ('........',
                 '........',
                 '...OX...',
                 '..XXXX..',
                 '..XXXOXX',
                 '...O.XO.',
                 '........',
                 '........')
        result = legal_moves(board, 'O')
        self.assertEqual({(2, 5), (2, 6), (3, 1), (3, 6), (4, 1), (5, 4), (6, 5)}, set(result))

    def test_finds_pass(self):
        board = ('........',
                 '........',
                 '........',
                 '...XOX..',
                 '........',
                 '........',
                 '........',
                 '........')
        result = legal_moves(board, 'X')
        self.assertEqual({'pass'}, set(result))

    def test_finds_no_legal_moves_before_board_full(self):
        board = ('OOOOOOOO',
                 'OOOOOOOO',
                 'OOOOOOOO',
                 'OOOOOOO.',
                 'OOOOOO..',
                 'OOOOOO.X',
                 'OOOOOOO.',
                 'OOOOOOOO')
        self.assertEqual(set(), set(legal_moves(board, 'X')))
        self.assertEqual(set(), set(legal_moves(board, 'O')))

    def test_finds_score(self):
        board = ('XXXXXXXX',
                 'XXXXXXXX',
                 'XXXXXXXX',
                 'XXXXXXXX',
                 'OOOOOOOX',
                 'OOOOOOOO',
                 'OOOOOOOO',
                 'OOOOOOOO')
        result = score(board)
        self.assertEqual(2, result)

    def test_finds_value_1(self):
        board = ('........',
                 '........',
                 '...OX...',
                 '..XXXX..',
                 '..XXXOXX',
                 '...O.XO.',
                 '........',
                 '........')
        v = value(board, 'O', 1)
        self.assertEqual(-2, v)

    def test_finds_value_2(self):
        board = ('........',
                 '........',
                 '...OX...',
                 '..XXXX..',
                 '..XXXOXX',
                 '...O.XO.',
                 '........',
                 '........')
        v = value(board, 'O', 2)
        self.assertEqual(7, v)

    def test_finds_best_move_shallow(self):
        board = ('XXXXXXXX',
                 'XXXXXXXX',
                 'XXXXXXXX',
                 '..OXXXXO',
                 'OOOOOOOO',
                 'OOOOOOOO',
                 'OOOOOOOO',
                 'XOOOOOOO')
        move = best_move(board, 'X', 1)
        self.assertEqual((3, 0), move)

    def test_finds_best_move_medium(self):
        board = ('........',
                 'OXOXX...',
                 '........',
                 '........',
                 '........',
                 '........',
                 '........',
                 'OX......')
        move = best_move(board, 'O', 2)
        self.assertEqual((7, 2), move)

    def test_finds_best_move_deep(self):
        board = ('........',
                 'OXOXX...',
                 '........',
                 '........',
                 '........',
                 '........',
                 '........',
                 'OX......')
        move = best_move(board, 'O', 3)
        self.assertEqual((1, 5), move)
