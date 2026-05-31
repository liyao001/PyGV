from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class GroupLabelConfig:
    start_track_idx: int
    end_track_idx: int
    label: str
    x: float = 0.02
    x_line_offset: float = 0.015

    def __post_init__(self):
        if self.end_track_idx <= self.start_track_idx:
            raise ValueError("end_track_idx must be >= start_track_idx")
        if self.start_track_idx < 0 or self.end_track_idx < 0:
            raise ValueError("start_track_idx must be >= 0")


@dataclass
class GroupLabels:
    label_configs: List[GroupLabelConfig] = field(default_factory=list)

    def add(self, label_config: GroupLabelConfig):
        self.label_configs.append(label_config)
