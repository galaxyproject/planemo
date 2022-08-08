import argparse
import ast
from typing import Optional

class DecoyParser(argparse.ArgumentParser):
    """
    Decoy class that is injected into code that initializes parser
    It removes dependency on the internal representation of actions
    """

    class Action:
        def __init__(self, scope, argument, action,
                     kwargs):
            self.argument = argument

            if action is None:
                action = "STORE"

            self.action = action.upper()
            self.kwargs = kwargs
            self.scope = scope

    class Section:
        def __init__(self, name=None, description=None):
            self.name = name
            self.description = description

    def __init__(self):
        self.tracked_actions = []
        self.default_scope = self.Section("default")
        super().__init__()

    def report_arguments_and_groups(self):
        return self.tracked_actions

    def save_action(self, scope, *args, **kwargs):
        self.tracked_actions.append(
            self.Action(scope, args[-1], kwargs.get("action",
                                                   "STORE"), kwargs))

    def add_argument(self, *args, **kwargs):
        self.save_action(self.default_scope, *args, **kwargs)
        return super().add_argument(*args, **kwargs)

    def add_argument_for_arg_group(self, scope, group):
        def add_argument(*args, **kwargs):
            self.save_action(scope, *args, **kwargs)
            return super(type(group), group).add_argument(*args, **kwargs)

        return add_argument

    def add_argument_group(self, *args, **kwargs):
        arg_group = super().add_argument_group(*args, **kwargs)
        if len(args) == 0:
            scope = self.Section(**kwargs)
        else:
            scope = self.Section(*args)

        arg_group.add_argument = self.add_argument_for_arg_group(scope,
                                                                 arg_group)
        return arg_group


def obtain_class_def() -> Optional[ast.ClassDef]:
    file = open(__file__, "r")
    module = ast.parse(file.read())
    file.close()
    return next((item for item in module.body if
                 type(item) is ast.ClassDef and item.name == "DecoyParser"),
                None)
