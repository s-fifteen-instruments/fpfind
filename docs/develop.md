# Development

## Local testing across different Python versions

Running tests with different Python can be easily done with `uv run`. For example, to test with Python 3.8 in an isolated environment:

```bash
uv run --python 3.8 --isolated --with-editable . --with pytest pytest
```

Editable installation is preferred to avoid issues from build caching.
The additional `pytest` dependency will need to be specified because `uv` does not recognize the `[tool.poetry.group]` table required by Poetry v1 for optional dependencies.

* Poetry v1 (`poetry-core~=1`) will reject unrecognized tables, including `[project.optional-dependencies]`
* Poetry v2 (`poetry-core~=2`) supports the `[project]` table but is available only for Python 3.9+

### Coverage and parallel test execution

Add `pytest-cov` for code coverage and `pytest-xdist` for parallel test execution mode:

```bash
uv run --python 3.8 --isolated --with-editable . --with pytest --with pytest-cov --with pytest-xdist pytest --cov=fpfind -n auto
```

## Lock files

Lock files can be generated in Poetry with `poetry lock` and uv with `uvx migrate-to-uv; uv lock`. For convenience, the latter is provided as a script in `scripts/generate-uv-lock.sh`.

These lock files however should not be committed to version control: differences in build will likely occur due to different build platforms, which can be a PITA to resolve (e.g. see this [SO answer](https://stackoverflow.com/a/61076546)). Instead, we build the library during CI and archive the generated lockfile as a build artifact (these can be found in the [individual CI workflow](https://github.com/s-fifteen-instruments/fpfind/actions)). This allows for reproducible builds and easier rollbacks when something breaks.
