import nox

@nox.session(python="3.11")
def tests(session):
    session.run("pip", "install", "pytest", external=True)
    session.run("pytest", "test_mesh_file.py", external=True)

@nox.session(python="3.11")
def docs(session):
    session.run("pip", "install", "sphinx", "myst_parser", "sphinx_autodoc_typehints", "sphinx_copybutton", external=True)
    session.run("make", "html", external=True)
