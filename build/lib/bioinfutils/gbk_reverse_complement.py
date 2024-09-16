from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

import argparse

# Function to reverse complement a GenBank file
def reverse_complement_genbank(input_file):
    # Parse the GenBank file
    records = SeqIO.parse(input_file, "genbank")

    # Create a list to store new features
    new_features = []

    # Iterate through each record in the GenBank file
    for record in records:
        # Reverse complement the sequence
        reversed_sequence = record.seq.reverse_complement()

        # Check if 'molecule_type' is present in annotations, set a default if not
        molecule_type = record.annotations.get("molecule_type", "DNA")

        # Create a new record with the reversed sequence
        reversed_record = record.__class__(
            id=record.id + "_reversed",
            name=record.name,
            description="Reversed complement of " + record.description,
            seq=reversed_sequence,
            features=[],
        )

        # Reverse complement features
        for feature in record.features:
            # Reverse complement gene features
            new_location = FeatureLocation(
                len(record.seq) - feature.location.end,
                len(record.seq) - feature.location.start,
                strand=-feature.location.strand,
            )
            new_feature = SeqFeature(
                location=new_location,
                type=feature.type,
                qualifiers=feature.qualifiers,
            )
            new_features.append(new_feature)

        # Add the new features to the reversed record
        reversed_record.features = new_features

        # Set the 'molecule_type' in the annotations of the reversed record
        reversed_record.annotations["molecule_type"] = molecule_type
        
        return reversed_record





def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', type=str, required=True)
    parser.add_argument('--output', '-o', type=str, required=True)
    
    args = parser.parse_args()

    revcomp_records = reverse_complement_genbank(args.input)
 
    SeqIO.write(revcomp_records, f"{args.output}", "genbank")
        

if __name__ == "__main__":
    main()
