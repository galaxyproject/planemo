"""
Module containing the parent class of Dicovery classes
"""
import ast
import abc
from typing import List, Tuple


class CustomVisitor(ast.NodeVisitor, abc.ABC):
    @abc.abstractmethod
    def report_findings(self) -> Tuple:
        pass

    def visit_and_report(self, source: ast.AST):
        self.visit(source)
        return self.report_findings()


class Discovery(CustomVisitor, abc.ABC):
    def __init__(self, actions: List[ast.AST]):
        self.actions = actions


def add_parents(tree: ast.AST):
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            child.parent = node


def create_module_tree_from_path(path: str) -> ast.Module:
    with open(path, mode="r", encoding="utf-8") as file:
        tree = ast.parse(file.read())
        add_parents(tree)
        return tree


def create_module_tree_from_str(text: str) -> ast.Module:
    tree = ast.parse(text)
    add_parents(tree)
    return tree
