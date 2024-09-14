# program_utils.py
import importlib

# def get_program(program_name):
#   """
#   Imports and returns the program module based on the given name.

#   Args:
#       program_name: The name of the program (e.g., 'program1').

#   Returns:
#       The imported program module or None if not found.
#   """
#   try:
#     program_module = importlib.import_module(f"api.{program_name}")
#     return program_module
#   except ModuleNotFoundError:
#     return None

def run_program(program_name):
  """
  Attempts to import the program and call its main function.

  Args:
      program_name: The name of the program (e.g., 'program1').
  """
  try:
    program_module = importlib.import_module(f"api.{program_name}")
    return program_module
  except ModuleNotFoundError:
    return None
  if program_module:
    program_module.main()
  else:
    print(f"Invalid program: {program_name}")