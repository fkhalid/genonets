
import os


def num_files_matches(ground_truth_dir, target_dir):
    num_gt_files = len(os.listdir(ground_truth_dir))
    num_target_files = len(os.listdir(target_dir))

    return num_gt_files == num_target_files
