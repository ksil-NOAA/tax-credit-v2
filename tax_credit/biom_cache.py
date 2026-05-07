#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Small LRU-style cache for repeated BIOM loads in mock evaluation (Phase 3).

Keys normalize paths with :func:`os.path.realpath` so the same file reached
via different symlinks shares one entry.
"""

from __future__ import annotations

from collections import OrderedDict
from os.path import realpath
from typing import Any, Callable, Hashable, Optional, TypeVar

T = TypeVar("T")


def normalized_table_path(table_fp: str) -> str:
    return realpath(table_fp)


def _taxa_cache_component(taxa_to_keep: Any) -> Any:
    if taxa_to_keep is None:
        return None
    if isinstance(taxa_to_keep, (list, tuple)):
        return tuple(taxa_to_keep)
    return (taxa_to_keep,)


def mount_observations_cache_key(
    table_fp: str,
    min_count: int,
    taxonomy_level: int,
    taxa_to_keep: Any,
    md_key: str,
    filter_obs: bool,
) -> tuple:
    """Cache key aligned with ``eval_framework.mount_observations`` call sites."""
    return (
        normalized_table_path(table_fp),
        int(min_count),
        int(taxonomy_level),
        _taxa_cache_component(taxa_to_keep),
        md_key,
        bool(filter_obs),
    )


def feature_table_cache_key(table_fp: str) -> tuple:
    return ("feature", normalized_table_path(table_fp))


class _NoBiomCache:
    """Passthrough: every lookup runs *factory* (no sharing)."""

    __slots__ = ()

    def get_or_put(self, key: Hashable, factory: Callable[[], T]) -> T:
        return factory()


class BiomTableCache:
    """OrderedDict-backed store with optional LRU eviction.

    Parameters
    ----------
    max_entries
        If ``None``, grow without eviction. If a positive int, drop the
        least-recently-used entry when a new key is inserted beyond this size.
    """

    __slots__ = ("_max", "_od")

    def __init__(self, max_entries: Optional[int] = None) -> None:
        if max_entries is not None and max_entries < 1:
            raise ValueError("max_entries must be None or >= 1")
        self._max = max_entries
        self._od: OrderedDict[Hashable, Any] = OrderedDict()

    def get_or_put(self, key: Hashable, factory: Callable[[], T]) -> T:
        if key in self._od:
            self._od.move_to_end(key)
            return self._od[key]
        val = factory()
        self._od[key] = val
        self._od.move_to_end(key)
        if self._max is not None and len(self._od) > self._max:
            self._od.popitem(last=False)
        return val

    def clear(self) -> None:
        self._od.clear()

    def __len__(self) -> int:
        return len(self._od)


NO_CACHE = _NoBiomCache()
