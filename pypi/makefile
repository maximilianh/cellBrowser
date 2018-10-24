all:
	rm -rf dist/*
	python2 setup.py sdist bdist_wheel
	python3 setup.py sdist bdist_wheel

up:
	twine upload dist/*

