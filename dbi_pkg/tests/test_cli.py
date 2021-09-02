import os
from click.testing import CliRunner
from dbinspector.cli import cli
from dbinspector.exceptions import QueryNotFoundError

TEST_FOLDER = os.path.dirname(__file__)


class TestCli:
    """Class for testing the CLI commands."""

    def test_cli(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli)
        assert result.exit_code == 0
        assert 'Entry method for the CLI.' in result.output  # Check if cli() called properly
        help_result = runner.invoke(cli, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output

    def test_parse(self):
        """Test the parse CLI command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['parse'])
        assert result.exit_code == 0  # cli command was successful
        help_result = runner.invoke(cli, ['parse', '--help'])
        assert "Parse the downloaded database data." in help_result.output

    def test_compare_successful(self):
        """Test the compare CLI command with a successful query."""
        runner = CliRunner()
        result = runner.invoke(cli, ['compare', '-q', 'Q9NY95'])
        assert result.exit_code == 0
        assert "symbol" in result.output

    def test_compare_outfile(self):
        """Test the compare CLI command with an output file given"""
        runner = CliRunner()
        outfile = "test-compare.tsv"
        result = runner.invoke(cli, ['compare', '-q', 'Q9NY95', '-o', outfile])
        assert os.path.isfile(outfile)
        os.remove(outfile)
        assert not os.path.isfile(outfile)

    def test_compare_unsuccessful(self):
        """Test the compare CLI command with an unsuccessful query."""
        runner = CliRunner()
        result = runner.invoke(cli, ['compare', '-q', 'blablabla'])
        assert isinstance(result.exception, QueryNotFoundError)
        assert result.exit_code == 1  # unsucessful

    def test_database_summary(self):
        """Test the database-summary CLI command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['database-summary'])
        assert result.exit_code == 0
        assert "Number of matches" in result.output

    def test_database_summary_outfile(self):
        """Test the database-summary CLI command with an output file given."""
        runner = CliRunner()
        outfile = "test-summary.tsv"
        result = runner.invoke(cli, ['database-summary', '-o', outfile])
        assert os.path.isfile(outfile)
        os.remove(outfile)
        assert not os.path.isfile(outfile)

    def test_check_age(self):
        """Test the check-age CLI command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['check-age', '-p', '-r'])
        assert result.exit_code == 0
        assert "Your raw refseq data is" in result.output and "Your raw uniprot data is" in result.output
        assert "Your parsed refseq data is" in result.output and "Your parsed uniprot data is" in result.output

    def test_clear_cache(self):
        """Test the clear-cache CLI command."""
        runner = CliRunner()
        result = runner.invoke(cli, ['clear-cache'], input="y")
        assert result.exit_code == 0
        assert "Clearing downloaded data files" in result.output

    def test_clear_cache_all(self):
        """Test the clear-cache CLI command with all param."""
        runner = CliRunner()
        result = runner.invoke(cli, ['clear-cache', '-a'], input="y")
        assert result.exit_code == 0
        assert "Clearing parsed data files" in result.output
        runner = CliRunner()
        # restore deleted data
        result2 = runner.invoke(cli, ['parse'])
        # sometimes output doesn't read the final output at the end, so check for first output message instead
        assert ("UniProt and RefSeq downloaded data parsed" in result2.output) or ("Downloading" in result2.output)
