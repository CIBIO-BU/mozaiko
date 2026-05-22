## mozaiko
VERSION=0.1.4

release:
	# Update version in pyproject.toml
	sed -i 's/^version *= *.*/version = "$(VERSION)"/' pyproject.toml
	sed -i 's/^__version__ =.*/__version__ ="$(VERSION)"/' src/mozaiko/mozaiko.py

	# Clean previous builds
	rm -rf dist build src/*.egg-info

	# Build package
	python3 -m build

	# Commit and tag
	git pull && \
	git add pyproject.toml && \
	git commit -e -m "New version $(VERSION)" . && \
	git tag -a "$(VERSION)" -m "v$(VERSION)" && \
	git push && git push --tags

	# Upload to PyPI (will prompt for API token)
	twine upload dist/*

setup:
	# run whole setup process
	bash conda_env_setup.sh

install-crabs:
	# run only the CRABS installation step
	bash -c "source conda_env_setup.sh && activate_env && install_crabs_release"