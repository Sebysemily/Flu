# Agent Temporary Scripts Rule

When creating temporary, exploratory, or one-off code that is NOT part of the main flow of the project's pipelines (e.g., scripts for testing, debugging, or fixing git conflicts):
- ALWAYS place these files in the `./scratch/` directory.
- DO NOT create temporary files in the root directory, `code/`, or any other main project directories.
- This ensures that temporary files are automatically ignored by git and keep the workspace clean.
