for fa in split_fa/*.fasta; do
  # extract first header without ">"
  header=$(head -n1 "$fa" | sed 's/^>//')
  # replace / and | with _
  clean_name=$(echo "$header" | tr '/|' '_')

  # output dir based on header
  out="my_outputs/$clean_name"
  mkdir -p "$out"

  # run inference (input stays as original file)
  python inference.py \
    --input_fas "$fa" \
    --device cuda:0 \
    --single_seq_pred True \
    --output_dir "$out" \
    --ckpt pretrained/RhoFold_pretrained.pt
done
