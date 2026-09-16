#!/usr/bin/env python

import unittest

from tax_credit.biom_cache import BiomTableCache, NO_CACHE, mount_observations_cache_key


class BiomTableCacheTests(unittest.TestCase):
    def test_no_cache_always_runs_factory(self):
        n = {'c': 0}

        def factory():
            n['c'] += 1
            return n['c']

        self.assertEqual(NO_CACHE.get_or_put('k', factory), 1)
        self.assertEqual(NO_CACHE.get_or_put('k', factory), 2)

    def test_cache_dedupes_same_key(self):
        c = BiomTableCache()
        n = {'c': 0}

        def factory():
            n['c'] += 1
            return 't'

        k = mount_observations_cache_key('/tmp/x', 0, 3, None, 'taxonomy', False)
        self.assertIs(c.get_or_put(k, factory), 't')
        self.assertIs(c.get_or_put(k, factory), 't')
        self.assertEqual(n['c'], 1)

    def test_lru_eviction(self):
        c = BiomTableCache(max_entries=2)
        c.get_or_put('a', lambda: 1)
        c.get_or_put('b', lambda: 2)
        c.get_or_put('a', lambda: 99)  # touch a
        c.get_or_put('c', lambda: 3)  # evicts b (oldest)
        self.assertEqual(c.get_or_put('a', lambda: 0), 1)
        self.assertEqual(c.get_or_put('c', lambda: 0), 3)
        self.assertEqual(len(c), 2)


if __name__ == '__main__':
    unittest.main()
