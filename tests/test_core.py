"""Tests for the core functionality."""

from datetime import datetime
import pytest
from pubmed_papers.core import PubMedFetcher, Author, Paper

def test_is_company_affiliation():
    """Test company affiliation detection."""
    fetcher = PubMedFetcher(email="test@example.com")
    
    # Test company affiliations
    assert fetcher._is_company_affiliation("Pfizer Inc.")
    assert fetcher._is_company_affiliation("Moderna Therapeutics")
    assert fetcher._is_company_affiliation("AstraZeneca Biotech")
    
    # Test academic affiliations
    assert not fetcher._is_company_affiliation("Harvard University")
    assert not fetcher._is_company_affiliation("Mayo Clinic")
    assert not fetcher._is_company_affiliation("National Institutes of Health")
    
    # Test edge cases
    assert not fetcher._is_company_affiliation("")
    assert not fetcher._is_company_affiliation(None)

def test_extract_email():
    """Test email extraction from text."""
    fetcher = PubMedFetcher(email="test@example.com")
    
    # Test valid emails
    assert fetcher._extract_email("Contact: test@example.com") == "test@example.com"
    assert fetcher._extract_email("Email: user.name@company.co.uk") == "user.name@company.co.uk"
    
    # Test no email
    assert fetcher._extract_email("No email here") is None
    assert fetcher._extract_email("") is None

def test_parse_publication_date():
    """Test publication date parsing."""
    fetcher = PubMedFetcher(email="test@example.com")
    
    # Test complete date
    assert fetcher._parse_publication_date({
        "Year": "2023",
        "Month": "12",
        "Day": "25"
    }) == datetime(2023, 12, 25)
    
    # Test partial date
    assert fetcher._parse_publication_date({
        "Year": "2023"
    }) == datetime(2023, 1, 1)
    
    # Test empty date
    assert fetcher._parse_publication_date({}) == datetime(1900, 1, 1)

def test_author_model():
    """Test Author model validation."""
    # Test valid author
    author = Author(
        name="John Doe",
        affiliation="Pfizer Inc.",
        email="john.doe@pfizer.com",
        is_corresponding=True,
        is_non_academic=True
    )
    assert author.name == "John Doe"
    assert author.affiliation == "Pfizer Inc."
    assert author.email == "john.doe@pfizer.com"
    assert author.is_corresponding
    assert author.is_non_academic
    
    # Test minimal author
    author = Author(name="Jane Smith")
    assert author.name == "Jane Smith"
    assert author.affiliation is None
    assert author.email is None
    assert not author.is_corresponding
    assert not author.is_non_academic

def test_paper_model():
    """Test Paper model validation."""
    author = Author(
        name="John Doe",
        affiliation="Pfizer Inc.",
        email="john.doe@pfizer.com",
        is_corresponding=True,
        is_non_academic=True
    )
    
    # Test valid paper
    paper = Paper(
        pubmed_id="12345",
        title="Test Paper",
        publication_date=datetime(2023, 1, 1),
        non_academic_authors=[author],
        company_affiliations=["Pfizer Inc."],
        corresponding_author_email="john.doe@pfizer.com"
    )
    assert paper.pubmed_id == "12345"
    assert paper.title == "Test Paper"
    assert paper.publication_date == datetime(2023, 1, 1)
    assert len(paper.non_academic_authors) == 1
    assert paper.non_academic_authors[0].name == "John Doe"
    assert paper.company_affiliations == ["Pfizer Inc."]
    assert paper.corresponding_author_email == "john.doe@pfizer.com" 
