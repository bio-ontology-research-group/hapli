from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

@dataclass
class AlignmentResult:
    """Stores the result of a single alignment."""
    feature_id: str
    feature_type: str
    mapq: int
    target_start: int
    target_end: int
    identity: float
    cigar: Optional[str] = None
    children: List['AlignmentResult'] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'feature_id': self.feature_id,
            'feature_type': self.feature_type,
            'mapq': self.mapq,
            'target_start': self.target_start,
            'target_end': self.target_end,
            'identity': self.identity,
            'cigar': self.cigar,
            'children': [child.to_dict() for child in self.children]
        }
