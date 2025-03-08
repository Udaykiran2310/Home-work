"""Tests for the command-line interface."""

from pathlib import Path
import pytest
from typer.testing import CliRunner
from pubmed_papers.cli import app

runner = CliRunner()

def test_help():
    """Test help output."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Options:" in result.stdout
    assert "--help" in result.stdout
    assert "--file" in result.stdout
    assert "--debug" in result.stdout

def test_missing_email():
    """Test that email is required."""
    result = runner.invoke(app, ["test query"])
    assert result.exit_code == 0  # Typer handles the prompt
    assert "email" in result.stdout.lower()

def test_debug_output(tmp_path):
    """Test debug output."""
    result = runner.invoke(
        app,
        [
            "test query",
            "--debug",
            "--email", "test@example.com",
            "--max-results", "1"
        ]
    )
    assert "Debug: Query: test query" in result.stdout

def test_output_file(tmp_path):
    """Test output to file."""
    output_file = tmp_path / "results.csv"
    result = runner.invoke(
        app,
        [
            "test query",
            "--file", str(output_file),
            "--email", "test@example.com",
            "--max-results", "1"
        ]
    )
    # Note: We can't test actual file contents without mocking PubMed API
    assert result.exit_code == 0 
