"""Core functionality for fetching and processing PubMed papers."""

from typing import List, Optional, Dict, Any
from datetime import datetime
import re
from Bio import Entrez

class PubMedFetcher:
    """Class to handle PubMed API interactions."""
    
    def __init__(self, email: str):
        """Initialize the PubMed fetcher with user email."""
        Entrez.email = email
        self.academic_keywords = {'university', 'college', 'institute', 'school',
                                'academia', 'laboratory', 'hospital', 'clinic',
                                'medical center', 'centre'}
        self.company_keywords = {'inc', 'corp', 'ltd', 'limited', 'llc',
                               'pharmaceutical', 'biotech', 'therapeutics',
                               'bioscience', 'laboratories'}
    
    def _is_company_affiliation(self, affiliation: str) -> bool:
        """Check if an affiliation is from a company."""
        if not affiliation:
            return False
        affiliation_lower = affiliation.lower()
        
        # Check for academic keywords
        if any(kw in affiliation_lower for kw in self.academic_keywords):
            return False
            
        # Check for company keywords
        return any(kw in affiliation_lower for kw in self.company_keywords)
    
    def _extract_email(self, text: str) -> Optional[str]:
        """Extract email address from text."""
        if not text:
            return None
        email_pattern = r'[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}'
        match = re.search(email_pattern, text)
        return match.group(0) if match else None
    
    def _parse_publication_date(self, date_dict: Dict[str, Any]) -> datetime:
        """Parse PubMed publication date into datetime object."""
        year = int(date_dict.get('Year', 1900))
        month = int(date_dict.get('Month', 1))
        day = int(date_dict.get('Day', 1))
        return datetime(year, month, day)
    
    def fetch_papers(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Fetch papers from PubMed based on the query.
        
        Args:
            query: PubMed search query
            max_results: Maximum number of results to return
            
        Returns:
            List of paper dictionaries
        """
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        papers = []
        
        # Fetch details for each paper
        for pubmed_id in record["IdList"]:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="xml")
            records = Entrez.parse(handle)
            
            for record in records:
                # Extract basic information
                title = record.get("ArticleTitle", "")
                pub_date = self._parse_publication_date(record.get("PubDate", {}))
                
                # Process authors and affiliations
                authors = []
                company_affiliations = set()
                corresponding_email = None
                
                for author in record.get("AuthorList", []):
                    name = f"{author.get('LastName', '')} {author.get('ForeName', '')}"
                    affiliation = author.get("Affiliation", "")
                    email = self._extract_email(affiliation)
                    
                    is_non_academic = self._is_company_affiliation(affiliation)
                    if is_non_academic:
                        company_affiliations.add(affiliation)
                    
                    author_data = {
                        "name": name.strip(),
                        "affiliation": affiliation,
                        "email": email,
                        "is_corresponding": False,
                        "is_non_academic": is_non_academic
                    }
                    
                    if email and not corresponding_email:
                        corresponding_email = email
                        author_data["is_corresponding"] = True
                    
                    if is_non_academic:
                        authors.append(author_data)
                
                if authors:  # Only include papers with non-academic authors
                    papers.append({
                        "pubmed_id": pubmed_id,
                        "title": title,
                        "publication_date": pub_date,
                        "non_academic_authors": authors,
                        "company_affiliations": list(company_affiliations),
                        "corresponding_author_email": corresponding_email
                    })
            
            handle.close()
        
        return papers 
