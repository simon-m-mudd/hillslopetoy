#!/usr/bin/env python

"""Tests for `hillslopetoy` package."""

import pytest

from click.testing import CliRunner

from hillslopetoy import hillslopetoy
from hillslopetoy import cli


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'hillslopetoy.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output

@pytest.fixture
def test_analytical():
    """
    This tests the anayltical solutions
    """
    
    x_half = hillslopetoy.set_profile_locations_half_length()
    x_full = hillslopetoy.set_profile_locations_constant()
    x_displace = hillslopetoy.displace_profile_locations_constant(x_half)
    
    z_half = hillslopetoy.ss_nonlinear_elevation(x_half)
    z_full = hillslopetoy.ss_nonlinear_elevation(x_full)
    z_displace = hillslopetoy.ss_nonlinear_elevation(x_displace)
    