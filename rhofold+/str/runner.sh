# Run RhoFold+ once per sequence
for fa in split_fa/*.fasta; do
  out="my_outputs/$(basename "${fa%.fasta}")"
  mkdir -p "$out"
  python inference.py \
    --input_fas "$fa" \
    --single_seq_pred True \
    --output_dir "$out" \
    --ckpt pretrained/RhoFold_pretrained.pt
done
