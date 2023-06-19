# Segmentation of cells using Cellpose 2.0

The processed differential phase contrast (DPC) images were segmented using [Cellpose 2.0](https://www.nature.com/articles/s41592-022-01663-4), where a pre-trained neural network model (Cyto) was improved using a human-in-the-loop approach. The human-in-the-loop approach used 125 image segments from multiple donors with user annotations for each image segment and a recursive training approach. The annotations and training were fine-tuned to segment red blood cells (healthy/round and sickle shapes) and to avoid segmentation of cells at the edges of the images and overlapping cells.

