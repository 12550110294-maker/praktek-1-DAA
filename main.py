# main.py
import csv
import os
from code05_proyek1.dataset.DNA_analyzer import DNAPattern, DNAAnalyzer


class Main:
    """Kelas utama pengesekusian program (Driver Class)."""
    
    @staticmethod
    def run():
        # Path file dataset CSV (naik 1 folder dari proyek 1daa)
        dataset_path = os.path.join(os.path.dirname(__file__), "..", "case1_dino_dna.csv")
        
        reference_dna = ""
        dna_samples = []

        # 1. Membaca File Dataset CSV
        if os.path.exists(dataset_path):
            with open(dataset_path, mode='r', encoding='utf-8') as file:
                reader = csv.reader(file)
                lines = [row for row in reader if row]
                
                if lines:
                    # Cek apakah baris pertama header atau data
                    start_idx = 0
                    if lines[0][0].strip().lower() in ["no", "id", "number"]:
                        start_idx = 1  # Lewati header
                    
                    # Baris data CSV (Format: No, Nama_Latin, Angka, DNA_Sequence)
                    for row in lines[start_idx:]:
                        if len(row) >= 4:
                            try:
                                num = int(row[0].strip())
                            except ValueError:
                                num = 0
                            species = row[1].strip()
                            pattern_str = row[3].strip()  # Kolom ke-4 adalah Sekuens DNA
                            
                            dna_samples.append(DNAPattern(pattern_str, species, num))
                    
                    # Jika ada data genom referensi terpisah, gunakan string referensi panjang
                    # Sebagai contoh dasar genom referensi untuk pengujian:
                    if dna_samples:
                        # Ambil salah satu pattern sebagai bagian dari genom referensi
                        reference_dna = dna_samples[0].pattern + "AGCTTCGGGCAACTTAACCTGGCACAAACCGTTCAGTCAG"

        # 2. Fallback jika file tidak terbaca
        if not reference_dna or not dna_samples:
            reference_dna = "AGCTTCGGGCAACTTAACCTGGCACAAACCGTTCAGTCAG"
            dna_samples = [
                DNAPattern("TTTTTTTTTTTTTTTT", "", 0),
                DNAPattern("CTGGCACAAACCGTTG", "Mapusaurus roseae", 90),
                DNAPattern("CTTCGGGCAACTTAAC", "Plateosaurus trossingensis", 45)
            ]

        # 3. Inisialisasi Analyzer dan Eksekusi KMP + Merge Sort
        analyzer = DNAAnalyzer(reference_dna)
        results = analyzer.process_patterns(dna_samples)

        # 4. Cetak Hasil Pemetaan (Tampilkan Top 5 & Bottom 3 agar tidak terlalu panjang)
        print("=== HASIL PEMETAAN FRAGMEN DNA (TERURUT SKOR TERTINGGI) ===")
        print(f"Total Fragmen Diperiksa: {len(results)} item\n")
        
        # Tampilkan 5 Teratas
        for idx, item in enumerate(results[:5], 1):
            print(f"{idx}. Pola: \"{item.pattern}\"")
            print(f"   Nama Latin: {item.species_name}")
            print(f"   Nomor: {item.number if item.number != 0 else ''}")
            print(f"   Indeks: {item.index if item.index != -1 else ''}")
            print(f"   Skor Kemiripan: {item.similarity:.1f}% ({item.matched_chars}/{len(item.pattern)})")
            print(f"   Status: {item.status}")
            print()

        if len(results) > 5:
            print(f"... [{len(results) - 5} data lainnya telah diurutkan] ...\n")

        print("Metrik Eksekusi:")
        print(f"Total Operasi Dasar (Perbandingan Karakter): {analyzer.matcher.operation_count} operasi")


if __name__ == "__main__":
    Main.run()