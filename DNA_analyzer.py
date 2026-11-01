# dna_analyzer.py

class DNAPattern:
    """Kelas untuk memodelkan data fragmen DNA dan hasil pencariannya."""
    def __init__(self, pattern: str, species_name: str = "", number: int = 0):
        self.pattern = pattern
        self.species_name = species_name
        self.number = number
        self.index = -1
        self.similarity = 0.0
        self.matched_chars = 0
        self.status = "Not Found"

    def __repr__(self):
        return f"DNAPattern({self.pattern}, Index={self.index}, Score={self.similarity:.1f}%, Status={self.status})"


class KMPMatcher:
    """Kelas pendeteksi string menggunakan Algoritma Knuth-Morris-Pratt dengan End-Mutation Tolerance."""
    def __init__(self):
        self.operation_count = 0  # Metrik menghitung perbandingan karakter

    def _compute_lps(self, pattern: str) -> list:
        """Membuat tabel Longest Prefix Suffix (LPS) untuk KMP."""
        m = len(pattern)
        lps = [0] * m
        length = 0
        i = 1

        while i < m:
            self.operation_count += 1
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    def search_exact(self, text: str, pattern: str) -> int:
        """Mencari kemunculan eksak (100% cocok) menggunakan KMP."""
        n = len(text)
        m = len(pattern)
        if m == 0 or n < m:
            return -1

        lps = self._compute_lps(pattern)
        i = 0
        j = 0

        while i < n:
            self.operation_count += 1
            if text[i] == pattern[j]:
                i += 1
                j += 1

            if j == m:
                return i - j  # Indeks awal ditemukan
            elif i < n and text[i] != pattern[j]:
                self.operation_count += 1
                if j != 0:
                    j = lps[j - 1]
                else:
                    i += 1
        return -1

    def search_end_mutated(self, text: str, pattern: str) -> int:
        """Mencari kemunculan dengan toleransi mutasi 1 basa di paling ujung kanan."""
        n = len(text)
        m = len(pattern)
        if m <= 1 or n < m:
            return -1

        prefix_pattern = pattern[:-1]
        m_prefix = len(prefix_pattern)
        lps = self._compute_lps(prefix_pattern)

        i = 0
        j = 0

        while i < n:
            self.operation_count += 1
            if text[i] == prefix_pattern[j]:
                i += 1
                j += 1

            if j == m_prefix:
                start_idx = i - j
                end_idx = start_idx + m - 1
                
                if end_idx < n:
                    self.operation_count += 1
                    return start_idx
                
                j = lps[j - 1]
            elif i < n and text[i] != prefix_pattern[j]:
                self.operation_count += 1
                if j != 0:
                    j = lps[j - 1]
                else:
                    i += 1
        return -1


class MergeSorter:
    """Kelas pengurutan manual Merge Sort."""
    @staticmethod
    def sort(arr: list) -> list:
        if len(arr) <= 1:
            return arr

        mid = len(arr) // 2
        left = MergeSorter.sort(arr[:mid])
        right = MergeSorter.sort(arr[mid:])

        return MergeSorter._merge(left, right)

    @staticmethod
    def _merge(left: list, right: list) -> list:
        result = []
        i = j = 0

        while i < len(left) and j < len(right):
            # Prioritas 1: Skor Kemiripan tertinggi
            if left[i].similarity > right[j].similarity:
                result.append(left[i])
                i += 1
            elif left[i].similarity < right[j].similarity:
                result.append(right[j])
                j += 1
            else:
                # Prioritas 2: Indeks terkecil
                if left[i].index <= right[j].index and left[i].index != -1:
                    result.append(left[i])
                    i += 1
                else:
                    result.append(right[j])
                    j += 1

        while i < len(left):
            result.append(left[i])
            i += 1
        while j < len(right):
            result.append(right[j])
            j += 1

        return result


class DNAAnalyzer:
    """Kelas pengelola analisis genomika."""
    def __init__(self, reference_genome: str):
        self.reference_genome = reference_genome
        self.matcher = KMPMatcher()

    def process_patterns(self, patterns: list) -> list:
        for item in patterns:
            k = len(item.pattern)
            
            # Step 1: Exact Match (100%)
            exact_idx = self.matcher.search_exact(self.reference_genome, item.pattern)
            if exact_idx != -1:
                item.index = exact_idx
                item.matched_chars = k
                item.similarity = 100.0
                item.status = "Exact Match"
                continue

            # Step 2: End-Mutated Match (Mutasi 1 basa di ujung)
            mutated_idx = self.matcher.search_end_mutated(self.reference_genome, item.pattern)
            if mutated_idx != -1:
                item.index = mutated_idx
                item.matched_chars = k - 1
                item.similarity = ((k - 1) / k) * 100.0
                item.status = "End-Mutated"
                continue

            # Step 3: Not Found
            item.index = -1
            item.matched_chars = 0
            item.similarity = 0.0
            item.status = "Not Found"

        # Step 4: Urutkan Hasil menggunakan Merge Sort
        sorted_patterns = MergeSorter.sort(patterns)
        return sorted_patterns