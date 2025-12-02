##############################################################################
# (c) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Crufty tool to detect global variables in Fortran source. This is a stop-gap
to aid the eradication of globals until Stylist can do it properly.
"""
from dataclasses import dataclass
from errno import ENOENT
from logging import getLogger
from os import strerror
from pathlib import Path
from re import compile as re_compile
from typing import Any, Callable, List, Optional, Sequence, Tuple

from fparser.common.readfortran import FortranFileReader  # type: ignore
from fparser.two.Fortran2003 import (
    Attr_Spec,  # type: ignore
    Declaration_Type_Spec,
    Entity_Decl,
    Function_Reference,
    Initialization,
    Intrinsic_Type_Spec,
    Main_Program,
    Module,
    Name,
    Program,
    Type_Declaration_Stmt,
)
from fparser.two.parser import ParserFactory  # type: ignore
from fparser.two.utils import get_child, walk, FortranSyntaxError  # type: ignore


# pylint: disable=too-few-public-methods
class Dirt:
    """
    Describes a dirty file.
    """

    def __init__(self, line_number: int, fortran_type: str, variable_name: str):
        self.line_number = line_number
        self.fortran_type = fortran_type.lower()
        self.variable_name = variable_name


# pylint: disable=too-few-public-methods
class DirtyFile:
    """
    List of dirty files.
    """

    # pylint: disable = redefined-outer-name
    def __init__(self, filename: Path) -> None:
        self.filename = filename
        self.dirt: List[Dirt] = []

    def __lt__(self, other: Any):
        """
        Compares this object against another for "less than" relationship.

        This is done based on filename.
        """
        if not isinstance(other, DirtyFile):
            message = (
                "Can only compare DirtyFile against other DirtyFile " "objects."
            )
            raise ValueError(message)
        return self.filename < other.filename

    def add_dirt(self, line_number: int, fortran_type: str, variable_name: str):
        """
        Extends the list of dirty files.
        """
        self.dirt.append(Dirt(line_number, fortran_type, variable_name))


@dataclass
class Entity:
    """
    Information about a single declared entity.
    """

    name: str
    initialised: Initialization


@dataclass
class Declaration:
    """
    Information about a declaration.
    """

    line_number: int
    fortran_type: str
    attributes: List[str]
    entities: List[Entity]


def __find_declarations(
    root: Program,
    dirty_file: DirtyFile,
    handlers: Sequence[Callable[[DirtyFile, Program, Declaration], None]],
) -> None:
    for declaration in walk(root, Type_Declaration_Stmt):
        if isinstance(declaration.parent.parent, Main_Program):
            continue

        entities = [
            Entity(
                str(get_child(entity, Name)), get_child(entity, Initialization)
            )
            for entity in walk(declaration, Entity_Decl)
        ]

        intrinsic_type = get_child(declaration, Intrinsic_Type_Spec)
        user_type = get_child(declaration, Declaration_Type_Spec)
        actual_type = intrinsic_type or get_child(user_type, Name)

        declaration_info = Declaration(
            line_number=declaration.item.span[0],
            fortran_type=str(actual_type),
            attributes=[
                str(attribute).lower()
                for attribute in walk(declaration, Attr_Spec)
            ],
            entities=entities,
        )

        for handler in handlers:
            handler(dirty_file, declaration.parent.parent, declaration_info)


def __process_file(
    filename: Path,
    tree_handler: Sequence[Callable[[DirtyFile, Program, Declaration], None]],
) -> Optional[DirtyFile]:
    file_tally = DirtyFile(filename)

    reader = FortranFileReader(str(filename))
    # pylint: disable = redefined-outer-name
    parser: Program = ParserFactory().create(std="f2008")
    tree = parser(reader)
    __find_declarations(tree, file_tally, tree_handler)

    # pylint: disable = no-else-return
    if len(file_tally.dirt) == 0:
        return None
    else:
        return file_tally


def __find_globals(
    dirty_file: DirtyFile, parent: Program, declaration: Declaration
) -> None:
    if not isinstance(parent, Module):
        return

    if "parameter" in declaration.attributes:
        return

    for entity in declaration.entities:
        dirty_file.dirt.append(
            Dirt(declaration.line_number, declaration.fortran_type, entity.name)
        )


def __find_explicit_saved(
    dirty_file: DirtyFile, parent: Program, declaration: Declaration
) -> None:
    if isinstance(parent, (Module, Main_Program)):
        return

    if "save" not in declaration.attributes:
        return

    for entity in declaration.entities:
        dirty_file.dirt.append(
            Dirt(declaration.line_number, declaration.fortran_type, entity.name)
        )


def __find_implicit_saved(
    dirty_file: DirtyFile, parent: Program, declaration: Declaration
) -> None:
    if isinstance(parent, (Module, Main_Program)):
        return

    # No need to check these as they are immutable and so intrinsically
    # global.
    #
    if "parameter" in declaration.attributes:
        return

    # No need to check these as they are explicitly saved
    #
    if "save" in declaration.attributes:
        return

    # todo: Currently we ignore pointer initialisation as we haven't decided
    #       what to do about it.
    #
    for entity in declaration.entities:
        if entity.initialised is None:
            continue

        if get_child(entity.initialised, Function_Reference) is None:
            dirty_file.dirt.append(
                Dirt(
                    declaration.line_number,
                    declaration.fortran_type,
                    entity.name,
                )
            )


__LOG_MESSAGE = "{filename}: {fortran_type}: {names}"
__FORTRAN_EXTENSION_PATTERN = re_compile(r"\.[FfXx]90")


# pylint: disable=redefined-outer-name
def entry(
    file_objects: List[Path],
) -> Tuple[List[DirtyFile], List[Path], List[Path]]:
    """
    Descend file tree processing files.
    """
    dirty_list: List[DirtyFile] = []  # pylint: disable=redefined-outer-name
    clean_list: List[Path] = []  # pylint: disable=redefined-outer-name
    not_considered: List[Path] = []  # pylint: disable=redefined-outer-name

    while len(file_objects) > 0:
        file_object = file_objects.pop()
        if not file_object.exists():
            raise FileNotFoundError(ENOENT, strerror(ENOENT), file_object)
        if file_object.is_dir():
            getLogger("occupyfortran").debug("Descending into %s", file_object)
            file_objects.extend(file_object.iterdir())
        else:  # Object is a file
            if __FORTRAN_EXTENSION_PATTERN.match(file_object.suffix):
                try:
                    getLogger("occupyfortran").debug("Processing %s", file_object)
                    report = __process_file(
                        file_object,
                        [
                            __find_globals,
                            __find_explicit_saved,
                            __find_implicit_saved,
                        ],
                    )
                except FortranSyntaxError:
                    getLogger("occupyfortran").warning(f"Syntax error raised from fparser in {file_object}.")

                if report is None:
                    clean_list.append(file_object)
                else:  # File has dirt
                    if report is not None:
                        dirty_list.append(report)
            else:
                getLogger("occupyfortran").debug("Ignoring %s", file_object)
                not_considered.append(file_object)

    return dirty_list, clean_list, not_considered
