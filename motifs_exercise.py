import itertools
import random
from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO


# Task 1: Warm-up functions
lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

#motifs = ["ACGT", "ATGT", "CCGA"]

def count_matrix(motifs):
    #{"A": [2, 0, 0, 1], "C": [1, 2, 0, 0], "G": [0, 0, 3, 0], "T": [0, 1, 0, 2]}
    l = len(motifs[0])
    counts = {"A": [0]*l,
              "T": [0]*l,
              "C": [0]*l,
              "G": [0]*l,}
    for motif in motifs:
        for column, base in enumerate(motif):
            counts[base][column] += 1
    return counts


def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    total_score = 0

    for column in range(l):
        total_score += max(counts[base][column] for base in "ACGT")
    return total_score


def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])

    result = []
    for column in range(l):
        best_base = "A"
        best_count = counts["A"][column]
        for base in ["C", "G", "T"]:
            if counts[base][column] > best_count:
                best_base = base
                best_count = counts[base][column]
        result.append(best_base)
    return "".join(result)


def hamming_distance(a, b):
    distance = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    return distance


def total_distance(pattern, sequences):
    l = len(pattern)
    total_distance = 0

    for sequence in sequences:
        min_dist_in_seq = float("inf")

        for i in range(len(sequence) - l + 1):
            window = sequence[i: i + l]
            dist = hamming_distance(pattern, window)

            if dist < min_dist_in_seq:
                min_dist_in_seq = dist

        total_distance += min_dist_in_seq
    return total_distance

if __name__ == "__main__":
    red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]
    best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]
    print(score(red))  # 26
    print(consensus(best), score(best))  # AGATAG 34
    print(hamming_distance("TAAGTT", "TGAATT"))  # 2
    print(total_distance("TGCGTT", lecture_dna))  # 13



# Task 2: MotifProfile class

class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])             # store motif length l
        t = len(motifs)                     # number of sequences
        counts = count_matrix(motifs)       # counts the count matrix from task 1
        denominator = t + 4 * pseudocount

        # build position probability matrix (PPM) with pseudocounts
        self.ppm = {
            base: [
                (counts[base][col] + pseudocount) / denominator
                for col in range(self.l)
            ]
            for base in "ACGT"
        }

    def lmer_probability(self, lmer):
        """Multiply probabilities across all positions in the l-mer"""
        prob = 1.0
        for col, base in enumerate(lmer):
            prob *= self.ppm[base][col]
        return prob

    def most_probable_lmer(self, sequence):
        """Keep the track of the best window and its probability"""
        best_prob = -1.0
        best_window = ""

        # a window of length l sliding across the sequence
        for i in range(len(sequence) - self.l + 1):
            window = sequence[i : i + self.l]
            prob = self.lmer_probability(window)

            if prob > best_prob:
                best_prob = prob
                best_window = window

        return best_window

    def consensus(self):
        """Pick the most probable base for each column"""
        result = []
        for col in range(self.l):
            best_base = "A"
            best_prob = self.ppm["A"][col]
            for base in ["C", "G", "T"]:
                if self.ppm[base][col] > best_prob:
                    best_base = base
                    best_prob = self.ppm[base][col]
            result.append(best_base)
        return "".join(result)


if __name__ == "__main__":
    profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
    print(profile.consensus())  # ATGCGTA
    print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

    two = MotifProfile(["GTAC", "TTAA"])
    print(two.most_probable_lmer("ACTGGATGACCC"))  # TGAC
    print(round(two.lmer_probability("TGAC"), 4))  # 0.0093


instances = [Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]]
bio = motifs.create(instances)
bio.pseudocounts = 1

print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]


# The random module
rng = random.Random(1)                 # a random number generator with seed 1
i = rng.randint(0, 50)                 # random integer, 0 <= i <= 50 (both ends included!)
lmer = rng.choice(["ACG", "CGT", "GTA"])   # one random item of a list

# Task 3: MotifFinder class

sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]

class MotifFinder:

    def __init__(self, sequences, l, seed=None):
        # Store sequences and motif length
        self.sequences = sequences
        self.l = l

        self.rng = random.Random(seed)

        self.windows = [
            [seq[i : i + l] for i in range(len(seq) - l + 1)]
            for seq in self.sequences
        ]

    def total_distance(self, pattern):
        """Sum of minimal Hamming distances from pattern to all sequences."""
        total_dist = 0
        for seq_windows in self.windows:
            # find the closest window in the precomputed list
            min_dist = min(
                hamming_distance(pattern, window) for window in seq_windows
            )
            total_dist += min_dist
        return total_dist

    def median_string(self):
        """Find the pattern of length l that minimizes total distance across all 4^l combinations."""
        best_pattern = ""
        min_total_dist = float("inf")

        # generate all 4^l combinations of bases
        for combo in itertools.product("ACGT", repeat=self.l):
            pattern = "".join(combo)
            dist = self.total_distance(pattern)

            # keep the pattern with the lowest total distance
            if dist < min_total_dist:
                min_total_dist = dist
                best_pattern = pattern

        return best_pattern, min_total_dist

    def randomized_search(self):
        """Perform a single run of the randomized motif search."""
        # 1. randomly pick one l-mer from each sequence
        current_motifs = [
            self.rng.choice(seq_windows) for seq_windows in self.windows
        ]
        current_score = score(current_motifs)

        while True:
            # 2. build profile from current motifs with pseudocount=1
            profile = MotifProfile(current_motifs, pseudocount=1)

            # 3. find the profile-most-probable l-mer in every sequence
            new_motifs = [profile.most_probable_lmer(seq) for seq in self.sequences]
            new_score = score(new_motifs)

            # 4. if new motifs have a strictly higher score, continue
            if new_score > current_score:
                current_motifs = new_motifs
                current_score = new_score
            else:
                # else, local optimum reached; stop and return
                return current_motifs, current_score

    def best_of(self, runs):
        """Run randomized_search `runs` times and return the best motifs and score."""
        best_motifs = None
        best_score = -1

        for _ in range(runs):
            motifs_found, run_score = self.randomized_search()
            if run_score > best_score:
                best_score = run_score
                best_motifs = motifs_found

        return best_motifs, best_score


if __name__ == "__main__":
    from Bio import SeqIO

    # 1. Kontrola na základních datech (l = 6)
    print("\n Base data")
    finder_lec = MotifFinder(lecture_dna, 6, seed=1)
    print("Median string:", finder_lec.median_string())  # ('AGATAG', 2)
    motifs_lec, sc_lec = finder_lec.best_of(100)
    print(
        "Randomized:", consensus(motifs_lec), sc_lec
    )  # AGATAG 34

    # 2. Kontrola na souboru planted_motif.fasta (l = 7)
    print("\n Planted_motif.fasta")
    seqs = [
        str(rec.seq) for rec in SeqIO.parse("planted_motif.fasta", "fasta")
    ]
    finder_planted = MotifFinder(seqs, 7, seed=1)
    print("Median string:", finder_planted.median_string())  # ('GCTAAAG', 10)
    motifs_p, sc_p = finder_planted.best_of(100)
    print("Randomized:", consensus(motifs_p), sc_p)  # GCTAAAG 60


