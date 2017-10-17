"""Setuptools installation script for statsnerds repo."""

from setuptools import setup, find_packages


setup(
	name='statsnerds-pubs-2017',
	version='0.1',
	description='Team Stats Nerds pipeline for PUBS 2017',
	packages=find_packages('statsnerds'),
	install_requires=[
		'Biopython>=1.69',
		'numpy>=1.12'
	],
	include_package_data=True,
)
