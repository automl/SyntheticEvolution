from pathlib import Path
import pandas as pd
import numpy as np
import ast

# neurips ids come from running python results_summary_neurips.py if no ids are there, uncomment ids below for first run, then results_summary_neurips.py, then rerun preprocessing with only usable ids...
try:
    from neurips_ids import *  # noqa: F403
except ImportError:
    pass

# '6TZ2',

# # # 377 protein rna
# # protein_rna_ids = ["1A9N", "1DFU", "1EC6", "1EIY", "1EUQ", "1EUY", "1EXD", "1FEU", "1GAX", "1GTR", "1GTS", "1IVS", "1J1U", "1JBT", "1KUQ", "1LNG", "1MMS", "1O0B", "1O0C", "1OOA", "1P6V", "1QA6", "1QRS", "1QRT", "1QRU", "1QTQ", "1RPU", "1S03", "1SDS", "1SJ3", "1SJ4", "1SJF", "1U0B", "1UN6", "1VBX", "1VBY", "1VBZ", "1VC0", "1VC6", "1VFG", "1ZBH", "1ZDK", "1ZJW", "2AKE", "2ANN", "2B2G", "2BQ5", "2DER", "2DET", "2DR2", "2DR5", "2DR7", "2DR8", "2DR9", "2DRA", "2DRB", "2DU4", "2DU5", "2DU6", "2EZ6", "2F8K", "2F8S", "2F8T", "2GIC", "2HVY", "2IY5", "2NUE", "2NUF", "2OIH", "2OJ3", "2OZB", "2PJP", "2PLY", "2R8S", "2RD2", "2RE8", "2UWM", "2V3C", "2W2H", "2ZH1", "2ZH2", "2ZH3", "2ZH4", "2ZH5", "2ZH6", "2ZH7", "2ZH8", "2ZH9", "2ZHA", "2ZHB", "2ZKO", "2ZNI", "2ZXU", "2ZZN", "3A2K", "3AL0", "3AMU", "3EPJ", "3FOZ", "3FTF", "3HTX", "3IVK", "3IWN", "3K0J", "3KTW", "3NDB", "3NVI", "3OUY", "3OV7", "3OVA", "3OVS", "3PU4", "3R9W", "3R9X", "3SIU", "3TUP", "3U4M", "3U56", "3UMY", "3V7E", "3VJR", "3WBM", "3WC2", "3WQY", "4AL5", "4AL6", "4C7O", "4FTS", "4GHL", "4IG8", "4J39", "4JGN", "4JK0", "4JNG", "4KJI", "4KR6", "4KR9", "4L8H", "4M2Z", "4M30", "4M4O", "4NGF", "4NH6", "4OOG", "4P3E", "4PDB", "4PMI", "4PRF", "4QG3", "4QI2", "4QIL", "4QOZ", "4QQB", "4RWN", "4RWO", "4RWP", "4S3N", "4TUW", "4TUX", "4WC2", "4WC3", "4WC4", "4WC5", "4WC6", "4WC7", "4X0A", "4X0B", "4X4O", "4XCO", "4YB1", "4YCO", "4YVI", "4YVJ", "4YVK", "4YYE", "4ZLD", "4ZT9", "5CCB", "5D6G", "5DAR", "5DE5", "5DV7", "5F5F", "5F5H", "5HC9", "5JC7", "5MSF", "5NPM", "5UD5", "5UDZ", "5V6X", "5VW1", "5W0O", "5W1H", "5WLH", "5WT3", "5Y7M", "5ZW4", "6AAX", "6B14", "6B3K", "6BJG", "6BJH", "6BJV", "6BJY", "6CF2", "6CMN", "6DB8", "6DB9", "6DU4", "6E9E", "6E9F", "6G1S", "6G1X", "6GKH", "6H61", "6H66", "6IV6", "6KL9", "6KLB", "6NMA", "6NMC", "6NMD", "6O5F", "6P7M", "6SX0", "6SX2", "6SY6", "6TQB", "6U8D", "6VRB", "6VRC", "6W62", "6XH1", "6XH2", "6XH3", "6XJQ", "6XJW", "6XJY", "6XJZ", "6ZDQ", "6ZLC", "7BKQ", "7CTT", "7D7V", "7DLZ", "7DMQ", "7DOK", "7DTE", "7DWH", "7ED5", "7K9M", "7KA0", "7KAB", "7KRP", "7MLX", "7N0B", "7N0C", "7NGA", "7NQ4", "7ORN", "7OS0", "7OZQ", "7P0V", "7PMQ", "7QR3", "7QR4", "7R97", "7RE2", "7SLP", "7SWF", "7U2A", "7UO4", "7UO7", "7UOB", "7UOE", "7VG2", "7VG3", "7VTN", "7W0C", "7W0D", "7W0E", "7WAH", "7WKP", "7X34", "7X7A", "7X7R", "7X8A", "7XC7", "7XSP", "7XSQ", "7XSR", "7XSS", "7XT4", "7XW2", "7Y80", "7Y81", "7Y82", "7Y84", "7Y85", "7Y8T", "7Y8Y", "7Y9X", "7Y9Y", "7YN9", "7YNA", "7YNB", "7YNC", "7YND", "7ZOL", "7ZOQ", "8AF0", "8ASW", "8AW3", "8AXF", "8BGU", "8C4U", "8CTI", "8D1V", "8D4B", "8D8N", "8D97", "8D9E", "8D9F", "8D9G", "8D9H", "8D9I", "8DFV", "8DG7", "8DGA", "8DVS", "8EEX", "8EEY", "8FZR", "8G6R", "8GNA", "8GS2", "8GU6", "8H0S", "8H69", "8HIO", "8HKE", "8HMY", "8IDF", "8ISS", "8IUO", "8JTR", "8OPS", "8OPT", "8OST", "8PTX", "8PTZ", "8S9T", "8S9U", "8S9V", "8S9X", "8SQ9", "8SSW", "8T5S", "8URB", "8WM4", "8WMC", "8WMI", "8WML", "8YE6", "8Z8X", "8Z90", "9CMP", "9CPO", "9DTT", "9ENB", "9ENC", "9ENE", "9ENF", "9IIY", "9IIZ", "9IJ2"]

# # rna_monomers_ids = ["1Q9A", "1Z43", "2A43", "2EET", "2EEV", "2G9C", "2GIS", "2XNW", "2XO0", "2XO1", "2YGH", "3A3A", "3DVZ", "3FO4", "3GAO", "3GOG", "3GX2", "3GX5", "3IQN", "3IQP", "3L0U", "3LA5", "3SD3", "483D", "4B5R", "4FE5", "4FEJ", "4FEN", "4FEO", "4FEP", "4L81", "4LVV", "4LVW", "4LVY", "4LVZ", "4LW0", "4LX5", "4LX6", "4TZX", "4TZY", "4XNR", "5FK2", "5FK3", "5FK5", "5FKE", "5FKF", "5M0H", "6D3P", "6JXM", "6UBU", "6UC7", "6UC8", "6VUI", "6WJS", "6Y0Y", "7D81", "7D82", "7E9I", "7EAF", "7EEM", "7EFG", "7EFH", "7EOG", "7EOI", "7EOK", "7EOL", "7EOM", "7EON", "7EOO", "7EOP", "7LYJ", "7MKT", "7MKY", "7QSH", "7SAM", "7SXP", "7TZT", "7URI", "7URM", "7VFT", "8HB8", "8KEB", "8KED", "8Q6N", "8TKK", "8U5J", "8UPT", "8UPY", "9C6K", "9CXF", "9DE8"]

# # all_data_ids = protein_rna_ids + rna_monomers_ids

# all_data_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '5M0H', '7URM', '7EEM',
#        '6E80', '8HB8', '7EOG', '6WJS', '7D7W', '3IQP', '4LVY', '6UC7',
#        '6UC8', '4LVZ', '3A3A', '4LX5', '7EOP', '6PQ7', '5DA6', '4CS1',
#        '5C5W', '4FEP', '8Q6N', '3IQN', '7EOO', '7EOK', '4XNR', '7URI',
#        '2YGH', '7SXP', '6D3P', '2EEV', '3D0M', '4TZX', '4LVW', '8UPT',
#        '7EFG', '4EN5', '2XNW', '1Q9A', '387D', '8KED', '2XO1', '7E9I',
#        '3DVZ', '4C40', '3GX2', '6TB7', '2ZY6', '7EFH', '2EET', '2A43',
#        '3R4F', '5FK3', '5FKE', '7EOI', '4TZY', '4B5R', '9CXF', '7EOL',
#        '6UBU', '7EON', '7LYJ', '3GCA', '4L81', '2G9C', '4LVV', '4FEJ',
#        '9DE8', '7QSH', '4LW0', '4LX6', '4AOB', '3GX5', '4FEN', '7EAF',
#        '6TFE', '5FKF', '7D81', '4FE5', '3SD3', '3Q50', '5FK5', '8UPY',
#        '8U5J', '7TZT', '7SAM', '4RBQ', '8KEB', '1Z43', '483D', '7VFT',
#        '8TKK', '9C6K', '3GOG', '3LA5', '6Y0Y', '3ND3', '6TF1', '2GIS',
#        '4P8Z', '4FEO', '7D82', '7MKT', '6TF3', '3FO4', '7KD1', '6VUI',
#        '3L0U', '7EOM', '3GAO', '5FK2', '9ENB', '5VW1', '8Z90', '1GTS',
#        '3OV7', '2UWM', '3NDB', '4J1G', '8D9H', '3OVA', '2ZH9', '6XH1',
#        '5NPM', '4AL5', '7OS0', '1EIY', '7WKP', '6E9E', '2PJP', '2HVY',
#        '4RWP', '4JGN', '4M4O', '7Y80', '5UD5', '4BHH', '1VBX',
#        '6SY6', '8E28', '7P0V', '5WLH', '6DTD', '2ZH2', '3FTF', '2RD2',
#        '1SJ4', '7BKQ', '7PMQ', '5HSW', '3U56', '7D7V', '2DR8', '7Y8Y',
#        '3AMU', '8HMY', '6NMD', '4KJI', '7W0C', '8D9G', '7UOB', '7OG0',
#        '7DTE', '8C4U', '3WC2', '2F8T', '1KUQ', '3V7E', '7XSP', '5DE5',
#        '4QG3', '1J1U', '8AXF', '2ZHA', '2DET', '6TQB', '8DFV', '2F8S',
#        '1ZJW', '8OMR', '8ISS', '2DU6', '8D9F', '4NGF', '2DR9', '6HTU',
#        '6P7M', '3IVK', '1VC6', '7Y85', '6XH2', '8AW3', '7W0E', '6U8D',
#        '2NUF', '6O5F', '4PMI', '4XCO', '3R9W', '7NQ4', '1LNG', '6IV6',
#        '2OZB', '7U2A', '7XSS', '2ZHB', '5WT3', '8URB', '6SX0', '2IY5',
#        '4YVK', '7X7A', '4TUX', '4P3E', '2ZH6', '7WAH', '7XT4', '8H69',
#        '8SQ9', '7YNB', '1QRS', '6G1S', '8DVS', '6BJV', '9CPO', '6B14',
#        '8FZR', '1ZBH', '8DGA', '7CTT', '8HIO', '7ZOQ', '4FTS', '6BJG',
#        '8H0S', '8JTR', '3EPJ', '4PRF', '4IG8', '1EUY', '7Y82', '7NIC',
#        '7UOE', '5DV7', '7X8A', '4WC4', '4KR6', '5DAR', '7R97', '7OZQ',
#        '8YE6', '1JBT', '3PU4', '1VFG', '9ENF', '4JK0', '7KRP', '7ED5',
#        '1S03', '8PTZ', '7DWH', '5V6X', '8D4B', '3IWN', '8R3L', '7N0B',
#        '4RWN', '5F5F', '7KAB', '1VC0', '4WC3', '7DOK', '2ZH8', '1QRU',
#        '8WMI', '7SWF', '5MSF', '1IVS', '4QIL', '8WMC', '2DR7', '6VRC',
#        '8D97', '1GAX', '3BOY', '8IUO', '2W2H', '7K9M', '9ENE', '1O0C',
#        '3K0J', '7Y81', '4OOG', '8S9V', '6NMA', '3UMY', '6XJZ', '3HTX',
#        '2OJ3', '6DB8', '8GU6', '6B3K', '7MLX', '8CTI', '4J39', '5ZW4',
#        '6DB9', '8D8N', '1VBZ', '4QI2', '2ZZN', '8AF0', '8G6R', '2ZNI',
#        '6ZDQ', '4QQB', '6XH3', '4YYE', '4WC7', '3OUY', '1GTR', '2V3C',
#        '2DRA', '3VJR', '2ZXU', '7YNA', '7ORN', '1QRT', '5D6G', '8PTX',
#        '6ZLC', '8SSW', '1EC6', '8E29', '7VG2', '7Y8T', '4L8H', '8T5S',
#        '7QR4', '6GKH', '9IIY', '1FEU', '3AL0', '6CMN', '4YCO', '6V5B',
#        '4WJ4', '5F5H', '6W62', '4X0B', '3SIU', '3NVI', '8E40', '8WML',
#        '1U0B', '4M2Z', '3OVS', '7XC7', '8GS2', '9IJ2', '4YB1', '7VG3',
#        '8S9X', '4JNG', '2OIH', '1A9N', '4GHL', '2DR5', '3R9X', '7UO7',
#        '6BJY', '4C7O', '6BJH', '9CMP', '4WC5', '2ZH3', '8EEX', '7XW2',
#        '7XSR', '4PDB', '4S3N', '8D9E', '6XJY', '1UN6', '1SJF', '4YVJ',
#        '2DR2', '7SLP', '3TUP', '6H61', '1EUQ', '7Y9Y', '7UO4', '2DRB',
#        '3FOZ', '2ZH5', '8OST', '7XSQ', '4M30', '3U4M', '7DMQ', '8ASW',
#        '6KL9', '6IV9', '2DVI', '7X34', '7ZOL', '5Y7M', '1SJ3', #  '5H0R',
#        '7KA0', '3WQY', '7NGA', '9IIZ', '4AL6', '3WBM', '2NUE', '7N0C',
#        '1RPU', '4ZLD', '2DU4', '6E9F', '7YND', '6VRB', '9ENC', '2ZH7',
#        '4TUW', '6NMC', '8BGU', '1DFU', '9DTT', '2F8K', '8WM4', '4NH6',
#        '2DU5', '7W0D', '6H66', '2PLY', '7VTN', '8GNA', '5JC7', '7DLZ',
#        '7M5O', '2EZ6', '7YN9', '5HC9', '6WW6', '8HKE', '8EEY', '4KR9',
#        '4X0A', '7Y9X', '1MMS', '2ZKO', '2DER', '1QTQ', '8Z8X', '7QR3',
#        '6SX2', '1VBY', '7Y84', '1ZDK', '8DG7', '5CCB', '4X4O', '6CF2',
#        '6KLB', '1OOA', '5UDZ', '8OPT', '2R8S', '2RE8', '1QA6', '6DU4',
#        '8D9I', '8S9T', '5W1H', '6G1X', '2ANN', '5W0O', '3A2K', '2AKE',
#        '4WC6', '4WC2', '4ZT9', '7X7R', '4QOZ', '7RE2', '2B2G', '8OPS',
#        '6XJQ', '3AM1', '2GIC', '6XJW', '4YVI', '1EXD', '8IDF', '1SDS',
#        '2ZH1', '8D1V', '2ZH4', '7YNC', '6AAX', '3KTW', '1P6V', '4RWO',
#        '2BQ5', '8S9U', '1O0B']

# #  '6TZ2',

# all_data_orphan_ids = ['5M0H', '7URM', '7EEM', '6E80', '8HB8', '6PQ7', '5DA6', '4CS1',
#        '5C5W', '7URI', '7SXP', '3D0M', '7EFG', '1Q9A', '387D', '7E9I',
#        '3DVZ', '4C40', '2ZY6', '7EFH', '2A43', '3GCA', '7QSH', '3Q50',
#        '8U5J', '4RBQ', '483D', '7VFT', '8TKK', '6Y0Y', '3ND3', '7MKT',
#        '6VUI', '8Z90', '3OV7', '2UWM', '4J1G', '8D9H', '3OVA', '2ZH9',
#        '6XH1', '4AL5', '7OS0', '7WKP', '2PJP', '4RWP', '4JGN',
#        '4M4O', '7Y80', '4BHH', '6SY6', '8E28', '7P0V', '5WLH', '6DTD',
#        '2ZH2', '3FTF', '7BKQ', '5HSW', '2DR8', '7Y8Y', '6NMD', '4KJI',
#        '7W0C', '8D9G', '7UOB', '7OG0', '7DTE', '8C4U', '2F8T', '7XSP',
#        '5DE5', '8AXF', '2ZHA', '6TQB', '2F8S', '8D9F', '4NGF', '2DR9',
#        '6HTU', '6P7M', '7Y85', '6XH2', '2NUF', '6O5F', '4PMI', '3R9W',
#        '7NQ4', '6IV6', '2OZB', '7U2A', '7XSS', '2ZHB', '8URB', '6SX0',
#        '7X7A', '4TUX', '2ZH6', '7WAH', '7XT4', '8H69', '8SQ9', '7YNB',
#        '6G1S', '8DVS', '6BJV', '9CPO', '6B14', '1ZBH', '8DGA', '7CTT',
#        '8HIO', '7ZOQ', '4FTS', '6BJG', '8H0S', '8JTR', '4IG8', '7Y82',
#        '7NIC', '7UOE', '5DV7', '7X8A', '4KR6', '7R97', '7OZQ', '1JBT',
#        '3PU4', '1VFG', '4JK0', '7KRP', '7ED5', '1S03', '8D4B', '8R3L',
#        '7N0B', '4RWN', '5F5F', '7DOK', '2ZH8', '8WMI', '7SWF', '5MSF',
#        '4QIL', '8WMC', '2DR7', '8D97', '3BOY', '8IUO', '2W2H', '7Y81',
#        '4OOG', '8S9V', '6NMA', '3HTX', '8GU6', '4J39', '5ZW4', '8D8N',
#        '4QI2', '8AF0', '8G6R', '4QQB', '6XH3', '3OUY', '2DRA', '3VJR',
#        '7YNA', '7ORN', '6ZLC', '8SSW', '1EC6', '8E29', '7VG2', '7Y8T',
#        '4L8H', '8T5S', '6GKH', '9IIY', '1FEU', '6CMN', '5F5H', '6W62',
#        '3SIU', '3NVI', '8E40', '8WML', '4M2Z', '3OVS', '7XC7', '8GS2',
#        '9IJ2', '7VG3', '8S9X', '4JNG', '1A9N', '4GHL', '2DR5', '7UO7',
#        '6BJY', '4C7O', '6BJH', '9CMP', '2ZH3', '8EEX', '7XW2', '7XSR',
#        '4PDB', '4S3N', '8D9E', '6H61', '7Y9Y', '7UO4', '2DRB', '2ZH5',
#        '7XSQ', '4M30', '6KL9', '2DVI', '7ZOL', '7NGA', '9IIZ',  # '5H0R',
#        '4AL6', '3WBM', '2NUE', '7N0C', '1RPU', '4ZLD', '7YND', '6VRB',
#        '2ZH7', '4TUW', '6NMC', '1DFU', '2F8K', '8WM4', '4NH6', '6H66',
#        '2PLY', '7VTN', '8GNA', '5JC7', '7M5O', '2EZ6', '7YN9', '6WW6',
#        '8HKE', '8EEY', '4KR9', '7Y9X', '2ZKO', '8Z8X', '6SX2', '7Y84',
#        '1ZDK', '8DG7', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ', '6DU4',
#        '8D9I', '8S9T', '5W1H', '6G1X', '2ANN', '5W0O', '7X7R', '4QOZ',
#        '7RE2', '2B2G', '2GIC', '8IDF', '1SDS', '2ZH1', '8D1V', '2ZH4',
#        '7YNC', '6AAX', '4RWO', '2BQ5', '8S9U']

# all_data_non_orphan_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '7EOG', '6WJS', '7D7W',
#        '3IQP', '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5', '7EOP',
#        '4FEP', '8Q6N', '3IQN', '7EOO', '7EOK', '4XNR', '2YGH', '6D3P',
#        '2EEV', '4TZX', '4LVW', '8UPT', '4EN5', '2XNW', '8KED', '2XO1',
#        '3GX2', '6TB7', '2EET', '3R4F', '5FK3', '5FKE', '7EOI', '4TZY',
#        '4B5R', '9CXF', '7EOL', '6UBU', '7EON', '7LYJ', '4L81', '2G9C',
#        '4LVV', '4FEJ', '9DE8', '4LW0', '4LX6', '4AOB', '3GX5', '4FEN',
#        '7EAF', '6TFE', '5FKF', '7D81', '4FE5', '3SD3', '5FK5', '8UPY',
#        '7TZT', '7SAM', '8KEB', '1Z43', '9C6K', '3GOG', '3LA5', '6TF1',
#        '2GIS', '4P8Z', '4FEO', '7D82', '6TF3', '3FO4', '7KD1', '3L0U',
#        '7EOM', '3GAO', '5FK2', '9ENB', '5VW1', '1GTS', '3NDB', '5NPM',
#        '1EIY', '6E9E', '2HVY', '5UD5', '1VBX', '2RD2', '1SJ4', '7PMQ',
#        '3U56', '7D7V', '3AMU', '8HMY', '3WC2', '1KUQ', '3V7E', '4QG3',
#        '1J1U', '2DET', '8DFV', '1ZJW', '8OMR', '8ISS', '2DU6', '3IVK',
#        '1VC6', '8AW3', '7W0E', '6U8D', '4XCO', '1LNG', '5WT3', '2IY5',
#        '4YVK', '4P3E', '1QRS', '8FZR', '3EPJ', '4PRF', '1EUY', '4WC4',
#        '5DAR', '8YE6', '9ENF', '8PTZ', '7DWH', '5V6X', '3IWN', '7KAB',
#        '1VC0', '4WC3', '1QRU', '1IVS', '6VRC', '1GAX', '7K9M', '9ENE',
#        '1O0C', '3K0J', '3UMY', '6XJZ', '2OJ3', '6DB8', '6B3K', '7MLX',
#        '8CTI', '6DB9', '1VBZ', '2ZZN', '2ZNI', '6ZDQ', '4YYE', '4WC7',
#        '1GTR', '2V3C', '2ZXU', '1QRT', '5D6G', '8PTX', '7QR4', '3AL0',
#        '4YCO', '6V5B', '4WJ4', '4X0B', '1U0B', '4YB1', '2OIH', '3R9X',
#        '4WC5', '6XJY', '1UN6', '1SJF', '4YVJ', '2DR2', '7SLP', '3TUP',
#        '1EUQ', '3FOZ', '8OST', '3U4M', '7DMQ', '8ASW', '6IV9', '7X34',
#        '5Y7M', '1SJ3', '7KA0', '3WQY', '2DU4', '6E9F', '9ENC', '8BGU',
#        '9DTT', '2DU5', '7W0D', '7DLZ', '5HC9', '4X0A', '1MMS', '2DER',
#        '1QTQ', '7QR3', '1VBY', '5CCB', '8OPT', '2R8S', '2RE8', '1QA6',
#        '3A2K', '2AKE', '4WC6', '4WC2', '4ZT9', '8OPS', '6XJQ', '3AM1',
#        '6XJW', '4YVI', '1EXD', '3KTW', '1P6V', '1O0B']

# # '6TZ2',

# all_data_b2021_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '5M0H', '6E80', '6WJS',
#        '7D7W', '3IQP', '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5',
#        '6PQ7', '5DA6', '4CS1', '5C5W', '4FEP', '3IQN', '4XNR', '2YGH',
#        '6D3P', '2EEV', '3D0M', '4TZX', '4LVW', '4EN5', '2XNW', '1Q9A',
#        '387D', '2XO1', '3DVZ', '4C40', '3GX2', '6TB7', '2ZY6', '2EET',
#        '2A43', '3R4F', '5FK3', '5FKE', '4TZY', '4B5R', '6UBU', '7LYJ',
#        '3GCA', '4L81', '2G9C', '4LVV', '4FEJ', '4LW0', '4LX6', '4AOB',
#        '3GX5', '4FEN', '7EAF', '6TFE', '5FKF', '7D81', '4FE5', '3SD3',
#        '3Q50', '5FK5', '4RBQ', '1Z43', '483D', '3GOG', '3LA5', '6Y0Y',
#        '3ND3', '6TF1', '2GIS', '4P8Z', '4FEO', '7D82', '6TF3', '3FO4',
#        '7KD1', '6VUI', '3L0U', '3GAO', '5FK2', '5VW1', '1GTS', '3OV7',
#        '2UWM', '3NDB', '4J1G', '3OVA', '2ZH9', '6XH1', '5NPM', '4AL5',
#        '1EIY', '6E9E', '2PJP', '2HVY', '4RWP', '4JGN', '4M4O',
#        '5UD5', '4BHH', '1VBX', '6SY6', '5WLH', '6DTD', '2ZH2', '3FTF',
#        '2RD2', '1SJ4', '5HSW', '3U56', '7D7V', '2DR8', '3AMU', '6NMD',
#        '4KJI', '3WC2', '2F8T', '1KUQ', '3V7E', '5DE5', '4QG3', '1J1U',
#        '2ZHA', '2DET', '6TQB', '2F8S', '1ZJW', '2DU6', '4NGF', '2DR9',
#        '6HTU', '6P7M', '3IVK', '1VC6', '6XH2', '6U8D', '2NUF', '6O5F',
#        '4PMI', '4XCO', '3R9W', '7NQ4', '1LNG', '6IV6', '2OZB', '2ZHB',
#        '5WT3', '6SX0', '2IY5', '4YVK', '4TUX', '4P3E', '2ZH6', '1QRS',
#        '6G1S', '6BJV', '6B14', '1ZBH', '7CTT', '4FTS', '6BJG', '3EPJ',
#        '4PRF', '4IG8', '1EUY', '5DV7', '4WC4', '4KR6', '5DAR', '1JBT',
#        '3PU4', '1VFG', '4JK0', '7KRP', '1S03', '5V6X', '3IWN', '7N0B',
#        '4RWN', '5F5F', '7KAB', '1VC0', '4WC3', '2ZH8', '1QRU', '5MSF',
#        '1IVS', '4QIL', '2DR7', '6VRC', '1GAX', '3BOY', '2W2H', '7K9M',
#        '1O0C', '3K0J', '4OOG', '6NMA', '3UMY', '3HTX', '2OJ3', '6DB8',
#        '6B3K', '7MLX', '4J39', '5ZW4', '6DB9', '1VBZ', '4QI2', '2ZZN',
#        '2ZNI', '6ZDQ', '4QQB', '6XH3', '4YYE', '4WC7', '3OUY', '1GTR',
#        '2V3C', '2DRA', '3VJR', '2ZXU', '1QRT', '5D6G', '6ZLC', '1EC6',
#        '4L8H', '6GKH', '1FEU', '3AL0', '6CMN', '4YCO', '6V5B', '4WJ4',
#        '5F5H', '6W62', '4X0B', '3SIU', '3NVI', '1U0B', '4M2Z', '3OVS',
#        '4YB1', '4JNG', '2OIH', '1A9N', '4GHL', '2DR5', '3R9X', '6BJY',
#        '4C7O', '6BJH', '4WC5', '2ZH3', '4PDB', '4S3N', '1UN6', '1SJF',
#        '4YVJ', '2DR2', '3TUP', '6H61', '1EUQ', '2DRB', '3FOZ', '2ZH5',
#        '4M30', '3U4M', '7DMQ', '6KL9', '6IV9', '2DVI', '5Y7M', '1SJ3',
#        '7KA0', '3WQY', '4AL6', '3WBM', '2NUE', '7N0C', '1RPU',  # '5H0R', 
#        '4ZLD', '2DU4', '6E9F', '6VRB', '2ZH7', '4TUW', '6NMC', '1DFU',
#        '2F8K', '4NH6', '2DU5', '6H66', '2PLY', '5JC7', '7M5O', '2EZ6',
#        '5HC9', '4KR9', '4X0A', '1MMS', '2ZKO', '2DER', '1QTQ', '6SX2',
#        '1VBY', '1ZDK', '5CCB', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ',
#        '2R8S', '2RE8', '1QA6', '6DU4', '5W1H', '6G1X', '2ANN', '5W0O',
#        '3A2K', '2AKE', '4WC6', '4WC2', '4ZT9', '4QOZ', '2B2G', '3AM1',
#        '2GIC', '4YVI', '1EXD', '1SDS', '2ZH1', '2ZH4', '6AAX', '3KTW',
#        '1P6V', '4RWO', '2BQ5', '1O0B']

# # '6TZ2',

# all_data_b2021_orphan_ids = ['5M0H', '6E80', '6PQ7', '5DA6', '4CS1', '5C5W', '3D0M', '1Q9A',
#        '387D', '3DVZ', '4C40', '2ZY6', '2A43', '3GCA', '3Q50', '4RBQ',
#        '483D', '6Y0Y', '3ND3', '6VUI', '3OV7', '2UWM', '4J1G', '3OVA',
#        '2ZH9', '6XH1', '4AL5', '2PJP', '4RWP', '4JGN', '4M4O',
#        '4BHH', '6SY6', '5WLH', '6DTD', '2ZH2', '3FTF', '5HSW', '2DR8',
#        '6NMD', '4KJI', '2F8T', '5DE5', '2ZHA', '6TQB', '2F8S', '4NGF',
#        '2DR9', '6HTU', '6P7M', '6XH2', '2NUF', '6O5F', '4PMI', '3R9W',
#        '7NQ4', '6IV6', '2OZB', '2ZHB', '6SX0', '4TUX', '2ZH6', '6G1S',
#        '6BJV', '6B14', '1ZBH', '7CTT', '4FTS', '6BJG', '4IG8', '5DV7',
#        '4KR6', '1JBT', '3PU4', '1VFG', '4JK0', '7KRP', '1S03', '7N0B',
#        '4RWN', '5F5F', '2ZH8', '5MSF', '4QIL', '2DR7', '3BOY', '2W2H',
#        '4OOG', '6NMA', '3HTX', '4J39', '5ZW4', '4QI2', '4QQB', '6XH3',
#        '3OUY', '2DRA', '3VJR', '6ZLC', '1EC6', '4L8H', '6GKH', '1FEU',
#        '6CMN', '5F5H', '6W62', '3SIU', '3NVI', '4M2Z', '3OVS', '4JNG',
#        '1A9N', '4GHL', '2DR5', '6BJY', '4C7O', '6BJH', '2ZH3', '4PDB',
#        '4S3N', '6H61', '2DRB', '2ZH5', '4M30', '6KL9', '2DVI',  #  '5H0R',
#        '4AL6', '3WBM', '2NUE', '7N0C', '1RPU', '4ZLD', '6VRB', '2ZH7',
#        '4TUW', '6NMC', '1DFU', '2F8K', '4NH6', '6H66', '2PLY', '5JC7',
#        '7M5O', '2EZ6', '4KR9', '2ZKO', '6SX2', '1ZDK', '4X4O', '6CF2',
#        '6KLB', '1OOA', '5UDZ', '6DU4', '5W1H', '6G1X', '2ANN', '5W0O',
#        '4QOZ', '2B2G', '2GIC', '1SDS', '2ZH1', '2ZH4', '6AAX', '4RWO',
#        '2BQ5']

# all_data_b2021_non_orphan_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '6WJS', '7D7W', '3IQP',
#        '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5', '4FEP', '3IQN',
#        '4XNR', '2YGH', '6D3P', '2EEV', '4TZX', '4LVW', '4EN5', '2XNW',
#        '2XO1', '3GX2', '6TB7', '2EET', '3R4F', '5FK3', '5FKE', '4TZY',
#        '4B5R', '6UBU', '7LYJ', '4L81', '2G9C', '4LVV', '4FEJ', '4LW0',
#        '4LX6', '4AOB', '3GX5', '4FEN', '7EAF', '6TFE', '5FKF', '7D81',
#        '4FE5', '3SD3', '5FK5', '1Z43', '3GOG', '3LA5', '6TF1', '2GIS',
#        '4P8Z', '4FEO', '7D82', '6TF3', '3FO4', '7KD1', '3L0U', '3GAO',
#        '5FK2', '5VW1', '1GTS', '3NDB', '5NPM', '1EIY', '6E9E', '2HVY',
#        '5UD5', '1VBX', '2RD2', '1SJ4', '3U56', '7D7V', '3AMU', '3WC2',
#        '1KUQ', '3V7E', '4QG3', '1J1U', '2DET', '1ZJW', '2DU6', '3IVK',
#        '1VC6', '6U8D', '4XCO', '1LNG', '5WT3', '2IY5', '4YVK', '4P3E',
#        '1QRS', '3EPJ', '4PRF', '1EUY', '4WC4', '5DAR', '5V6X', '3IWN',
#        '7KAB', '1VC0', '4WC3', '1QRU', '1IVS', '6VRC', '1GAX', '7K9M',
#        '1O0C', '3K0J', '3UMY', '2OJ3', '6DB8', '6B3K', '7MLX', '6DB9',
#        '1VBZ', '2ZZN', '2ZNI', '6ZDQ', '4YYE', '4WC7', '1GTR', '2V3C',
#        '2ZXU', '1QRT', '5D6G', '3AL0', '4YCO', '6V5B', '4WJ4', '4X0B',
#        '1U0B', '4YB1', '2OIH', '3R9X', '4WC5', '1UN6', '1SJF', '4YVJ',
#        '2DR2', '3TUP', '1EUQ', '3FOZ', '3U4M', '7DMQ', '6IV9', '5Y7M',
#        '1SJ3', '7KA0', '3WQY', '2DU4', '6E9F', '2DU5', '5HC9', '4X0A',
#        '1MMS', '2DER', '1QTQ', '1VBY', '5CCB', '2R8S', '2RE8', '1QA6',
#        '3A2K', '2AKE', '4WC6', '4WC2', '4ZT9', '3AM1', '4YVI', '1EXD',
#        '3KTW', '1P6V', '1O0B']

# all_data_a2021_ids = ['7URM', '7EEM', '8HB8', '7EOG', '7EOP', '8Q6N', '7EOO', '7EOK',
#        '7URI', '7SXP', '8UPT', '7EFG', '8KED', '7E9I', '7EFH', '7EOI',
#        '9CXF', '7EOL', '7EON', '9DE8', '7QSH', '8UPY', '8U5J', '7TZT',
#        '7SAM', '8KEB', '7VFT', '8TKK', '9C6K', '7MKT', '7EOM', '9ENB',
#        '8Z90', '8D9H', '7OS0', '7WKP', '7Y80', '8E28', '7P0V', '7BKQ',
#        '7PMQ', '7Y8Y', '8HMY', '7W0C', '8D9G', '7UOB', '7OG0', '7DTE',
#        '8C4U', '7XSP', '8AXF', '8DFV', '8OMR', '8ISS', '8D9F', '7Y85',
#        '8AW3', '7W0E', '7U2A', '7XSS', '8URB', '7X7A', '7WAH', '7XT4',
#        '8H69', '8SQ9', '7YNB', '8DVS', '9CPO', '8FZR', '8DGA', '8HIO',
#        '7ZOQ', '8H0S', '8JTR', '7Y82', '7NIC', '7UOE', '7X8A', '7R97',
#        '7OZQ', '8YE6', '9ENF', '7ED5', '8PTZ', '7DWH', '8D4B', '8R3L',
#        '7DOK', '8WMI', '7SWF', '8WMC', '8D97', '8IUO', '9ENE', '7Y81',
#        '8S9V', '6XJZ', '8GU6', '8CTI', '8D8N', '8AF0', '8G6R', '7YNA',
#        '7ORN', '8PTX', '8SSW', '8E29', '7VG2', '7Y8T', '8T5S', '7QR4',
#        '9IIY', '8E40', '8WML', '7XC7', '8GS2', '9IJ2', '7VG3', '8S9X',
#        '7UO7', '9CMP', '8EEX', '7XW2', '7XSR', '8D9E', '6XJY', '7SLP',
#        '7Y9Y', '7UO4', '8OST', '7XSQ', '8ASW', '7X34', '7ZOL', '7NGA',
#        '9IIZ', '7YND', '9ENC', '8BGU', '9DTT', '8WM4', '7W0D', '7VTN',
#        '8GNA', '7DLZ', '7YN9', '6WW6', '8HKE', '8EEY', '7Y9X', '8Z8X',
#        '7QR3', '7Y84', '8DG7', '8OPT', '8D9I', '8S9T', '7X7R', '7RE2',
#        '8OPS', '6XJQ', '6XJW', '8IDF', '8D1V', '7YNC', '8S9U']

# all_data_a2021_orphan_ids = ['7URM', '7EEM', '8HB8', '7URI', '7SXP', '7EFG', '7E9I', '7EFH',
#        '7QSH', '8U5J', '7VFT', '8TKK', '7MKT', '8Z90', '8D9H', '7OS0',
#        '7WKP', '7Y80', '8E28', '7P0V', '7BKQ', '7Y8Y', '7W0C', '8D9G',
#        '7UOB', '7OG0', '7DTE', '8C4U', '7XSP', '8AXF', '8D9F', '7Y85',
#        '7U2A', '7XSS', '8URB', '7X7A', '7WAH', '7XT4', '8H69', '8SQ9',
#        '7YNB', '8DVS', '9CPO', '8DGA', '8HIO', '7ZOQ', '8H0S', '8JTR',
#        '7Y82', '7NIC', '7UOE', '7X8A', '7R97', '7OZQ', '7ED5', '8D4B',
#        '8R3L', '7DOK', '8WMI', '7SWF', '8WMC', '8D97', '8IUO', '7Y81',
#        '8S9V', '8GU6', '8D8N', '8AF0', '8G6R', '7YNA', '7ORN', '8SSW',
#        '8E29', '7VG2', '7Y8T', '8T5S', '9IIY', '8E40', '8WML', '7XC7',
#        '8GS2', '9IJ2', '7VG3', '8S9X', '7UO7', '9CMP', '8EEX', '7XW2',
#        '7XSR', '8D9E', '7Y9Y', '7UO4', '7XSQ', '7ZOL', '7NGA', '9IIZ',
#        '7YND', '8WM4', '7VTN', '8GNA', '7YN9', '6WW6', '8HKE', '8EEY',
#        '7Y9X', '8Z8X', '7Y84', '8DG7', '8D9I', '8S9T', '7X7R', '7RE2',
#        '8IDF', '8D1V', '7YNC', '8S9U']

# all_data_a2021_non_orphan_ids = ['7EOG', '7EOP', '8Q6N', '7EOO', '7EOK', '8UPT', '8KED', '7EOI',
#        '9CXF', '7EOL', '7EON', '9DE8', '8UPY', '7TZT', '7SAM', '8KEB',
#        '9C6K', '7EOM', '9ENB', '7PMQ', '8HMY', '8DFV', '8OMR', '8ISS',
#        '8AW3', '7W0E', '8FZR', '8YE6', '9ENF', '8PTZ', '7DWH', '9ENE',
#        '6XJZ', '8CTI', '8PTX', '7QR4', '6XJY', '7SLP', '8OST', '8ASW',
#        '7X34', '9ENC', '8BGU', '9DTT', '7W0D', '7DLZ', '7QR3', '8OPT',
#        '8OPS', '6XJQ', '6XJW']

# rna_monomers_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '5M0H', '7URM', '7EEM',
#        '6E80', '8HB8', '7EOG', '6WJS', '7D7W', '3IQP', '4LVY', '6UC7',
#        '6UC8', '4LVZ', '3A3A', '4LX5', '7EOP', '6PQ7', '5DA6', '4CS1',
#        '5C5W', '4FEP', '8Q6N', '3IQN', '7EOO', '7EOK', '4XNR', '7URI',
#        '2YGH', '7SXP', '6D3P', '2EEV', '3D0M', '4TZX', '4LVW', '8UPT',
#        '7EFG', '4EN5', '2XNW', '1Q9A', '387D', '8KED', '2XO1', '7E9I',
#        '3DVZ', '4C40', '3GX2', '6TB7', '2ZY6', '7EFH', '2EET', '2A43',
#        '3R4F', '5FK3', '5FKE', '7EOI', '4TZY', '4B5R', '9CXF', '7EOL',
#        '6UBU', '7EON', '7LYJ', '3GCA', '4L81', '2G9C', '4LVV', '4FEJ',
#        '9DE8', '7QSH', '4LW0', '4LX6', '4AOB', '3GX5', '4FEN', '7EAF',
#        '6TFE', '5FKF', '7D81', '4FE5', '3SD3', '3Q50', '5FK5', '8UPY',
#        '8U5J', '7TZT', '7SAM', '4RBQ', '8KEB', '1Z43', '483D', '7VFT',
#        '8TKK', '9C6K', '3GOG', '3LA5', '6Y0Y', '3ND3', '6TF1', '2GIS',
#        '4P8Z', '4FEO', '7D82', '7MKT', '6TF3', '3FO4', '7KD1', '6VUI',
#        '3L0U', '7EOM', '3GAO', '5FK2']

# rna_monomers_orphan_ids = ['5M0H', '7URM', '7EEM', '6E80', '8HB8', '6PQ7', '5DA6', '4CS1',
#        '5C5W', '7URI', '7SXP', '3D0M', '7EFG', '1Q9A', '387D', '7E9I',
#        '3DVZ', '4C40', '2ZY6', '7EFH', '2A43', '3GCA', '7QSH', '3Q50',
#        '8U5J', '4RBQ', '483D', '7VFT', '8TKK', '6Y0Y', '3ND3', '7MKT',
#        '6VUI']

# rna_monomers_non_orphan_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '7EOG', '6WJS', '7D7W',
#        '3IQP', '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5', '7EOP',
#        '4FEP', '8Q6N', '3IQN', '7EOO', '7EOK', '4XNR', '2YGH', '6D3P',
#        '2EEV', '4TZX', '4LVW', '8UPT', '4EN5', '2XNW', '8KED', '2XO1',
#        '3GX2', '6TB7', '2EET', '3R4F', '5FK3', '5FKE', '7EOI', '4TZY',
#        '4B5R', '9CXF', '7EOL', '6UBU', '7EON', '7LYJ', '4L81', '2G9C',
#        '4LVV', '4FEJ', '9DE8', '4LW0', '4LX6', '4AOB', '3GX5', '4FEN',
#        '7EAF', '6TFE', '5FKF', '7D81', '4FE5', '3SD3', '5FK5', '8UPY',
#        '7TZT', '7SAM', '8KEB', '1Z43', '9C6K', '3GOG', '3LA5', '6TF1',
#        '2GIS', '4P8Z', '4FEO', '7D82', '6TF3', '3FO4', '7KD1', '3L0U',
#        '7EOM', '3GAO', '5FK2']

# rna_monomers_b2021_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '5M0H', '6E80', '6WJS',
#        '7D7W', '3IQP', '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5',
#        '6PQ7', '5DA6', '4CS1', '5C5W', '4FEP', '3IQN', '4XNR', '2YGH',
#        '6D3P', '2EEV', '3D0M', '4TZX', '4LVW', '4EN5', '2XNW', '1Q9A',
#        '387D', '2XO1', '3DVZ', '4C40', '3GX2', '6TB7', '2ZY6', '2EET',
#        '2A43', '3R4F', '5FK3', '5FKE', '4TZY', '4B5R', '6UBU', '7LYJ',
#        '3GCA', '4L81', '2G9C', '4LVV', '4FEJ', '4LW0', '4LX6', '4AOB',
#        '3GX5', '4FEN', '7EAF', '6TFE', '5FKF', '7D81', '4FE5', '3SD3',
#        '3Q50', '5FK5', '4RBQ', '1Z43', '483D', '3GOG', '3LA5', '6Y0Y',
#        '3ND3', '6TF1', '2GIS', '4P8Z', '4FEO', '7D82', '6TF3', '3FO4',
#        '7KD1', '6VUI', '3L0U', '3GAO', '5FK2']

# rna_monomers_b2021_orphan_ids = ['5M0H', '6E80', '6PQ7', '5DA6', '4CS1', '5C5W', '3D0M', '1Q9A',
#        '387D', '3DVZ', '4C40', '2ZY6', '2A43', '3GCA', '3Q50', '4RBQ',
#        '483D', '6Y0Y', '3ND3', '6VUI']

# rna_monomers_b2021_non_orphan_ids = ['6UES', '6JXM', '2XO0', '7MKY', '6UET', '6WJS', '7D7W', '3IQP',
#        '4LVY', '6UC7', '6UC8', '4LVZ', '3A3A', '4LX5', '4FEP', '3IQN',
#        '4XNR', '2YGH', '6D3P', '2EEV', '4TZX', '4LVW', '4EN5', '2XNW',
#        '2XO1', '3GX2', '6TB7', '2EET', '3R4F', '5FK3', '5FKE', '4TZY',
#        '4B5R', '6UBU', '7LYJ', '4L81', '2G9C', '4LVV', '4FEJ', '4LW0',
#        '4LX6', '4AOB', '3GX5', '4FEN', '7EAF', '6TFE', '5FKF', '7D81',
#        '4FE5', '3SD3', '5FK5', '1Z43', '3GOG', '3LA5', '6TF1', '2GIS',
#        '4P8Z', '4FEO', '7D82', '6TF3', '3FO4', '7KD1', '3L0U', '3GAO',
#        '5FK2']

# rna_monomers_a2021_ids = ['7URM', '7EEM', '8HB8', '7EOG', '7EOP', '8Q6N', '7EOO', '7EOK',
#        '7URI', '7SXP', '8UPT', '7EFG', '8KED', '7E9I', '7EFH', '7EOI',
#        '9CXF', '7EOL', '7EON', '9DE8', '7QSH', '8UPY', '8U5J', '7TZT',
#        '7SAM', '8KEB', '7VFT', '8TKK', '9C6K', '7MKT', '7EOM']

# rna_monomers_a2021_orphan_ids = ['7URM', '7EEM', '8HB8', '7URI', '7SXP', '7EFG', '7E9I', '7EFH',
#        '7QSH', '8U5J', '7VFT', '8TKK', '7MKT']

# rna_monomers_a2021_non_orphan_ids = ['7EOG', '7EOP', '8Q6N', '7EOO', '7EOK', '8UPT', '8KED', '7EOI',
#        '9CXF', '7EOL', '7EON', '9DE8', '8UPY', '7TZT', '7SAM', '8KEB',
#        '9C6K', '7EOM']

# # '6TZ2',

# protein_rna_ids = ['9ENB', '5VW1', '8Z90', '1GTS', '3OV7', '2UWM', '3NDB', '4J1G',
#        '8D9H', '3OVA', '2ZH9', '6XH1', '5NPM', '4AL5', '7OS0', '1EIY',
#        '7WKP', '6E9E', '2PJP', '2HVY', '4RWP', '4JGN', '4M4O',
#        '7Y80', '5UD5', '4BHH', '1VBX', '6SY6', '8E28', '7P0V', '5WLH',
#        '6DTD', '2ZH2', '3FTF', '2RD2', '1SJ4', '7BKQ', '7PMQ', '5HSW',
#        '3U56', '7D7V', '2DR8', '7Y8Y', '3AMU', '8HMY', '6NMD', '4KJI',
#        '7W0C', '8D9G', '7UOB', '7OG0', '7DTE', '8C4U', '3WC2', '2F8T',
#        '1KUQ', '3V7E', '7XSP', '5DE5', '4QG3', '1J1U', '8AXF', '2ZHA',
#        '2DET', '6TQB', '8DFV', '2F8S', '1ZJW', '8OMR', '8ISS', '2DU6',
#        '8D9F', '4NGF', '2DR9', '6HTU', '6P7M', '3IVK', '1VC6', '7Y85',
#        '6XH2', '8AW3', '7W0E', '6U8D', '2NUF', '6O5F', '4PMI', '4XCO',
#        '3R9W', '7NQ4', '1LNG', '6IV6', '2OZB', '7U2A', '7XSS', '2ZHB',
#        '5WT3', '8URB', '6SX0', '2IY5', '4YVK', '7X7A', '4TUX', '4P3E',
#        '2ZH6', '7WAH', '7XT4', '8H69', '8SQ9', '7YNB', '1QRS', '6G1S',
#        '8DVS', '6BJV', '9CPO', '6B14', '8FZR', '1ZBH', '8DGA', '7CTT',
#        '8HIO', '7ZOQ', '4FTS', '6BJG', '8H0S', '8JTR', '3EPJ', '4PRF',
#        '4IG8', '1EUY', '7Y82', '7NIC', '7UOE', '5DV7', '7X8A', '4WC4',
#        '4KR6', '5DAR', '7R97', '7OZQ', '8YE6', '1JBT', '3PU4', '1VFG',
#        '9ENF', '4JK0', '7KRP', '7ED5', '1S03', '8PTZ', '7DWH', '5V6X',
#        '8D4B', '3IWN', '8R3L', '7N0B', '4RWN', '5F5F', '7KAB', '1VC0',
#        '4WC3', '7DOK', '2ZH8', '1QRU', '8WMI', '7SWF', '5MSF', '1IVS',
#        '4QIL', '8WMC', '2DR7', '6VRC', '8D97', '1GAX', '3BOY', '8IUO',
#        '2W2H', '7K9M', '9ENE', '1O0C', '3K0J', '7Y81', '4OOG', '8S9V',
#        '6NMA', '3UMY', '6XJZ', '3HTX', '2OJ3', '6DB8', '8GU6', '6B3K',
#        '7MLX', '8CTI', '4J39', '5ZW4', '6DB9', '8D8N', '1VBZ', '4QI2',
#        '2ZZN', '8AF0', '8G6R', '2ZNI', '6ZDQ', '4QQB', '6XH3', '4YYE',
#        '4WC7', '3OUY', '1GTR', '2V3C', '2DRA', '3VJR', '2ZXU', '7YNA',
#        '7ORN', '1QRT', '5D6G', '8PTX', '6ZLC', '8SSW', '1EC6', '8E29',
#        '7VG2', '7Y8T', '4L8H', '8T5S', '7QR4', '6GKH', '9IIY', '1FEU',
#        '3AL0', '6CMN', '4YCO', '6V5B', '4WJ4', '5F5H', '6W62', '4X0B',
#        '3SIU', '3NVI', '8E40', '8WML', '1U0B', '4M2Z', '3OVS', '7XC7',
#        '8GS2', '9IJ2', '4YB1', '7VG3', '8S9X', '4JNG', '2OIH', '1A9N',
#        '4GHL', '2DR5', '3R9X', '7UO7', '6BJY', '4C7O', '6BJH', '9CMP',
#        '4WC5', '2ZH3', '8EEX', '7XW2', '7XSR', '4PDB', '4S3N', '8D9E',
#        '6XJY', '1UN6', '1SJF', '4YVJ', '2DR2', '7SLP', '3TUP', '6H61',
#        '1EUQ', '7Y9Y', '7UO4', '2DRB', '3FOZ', '2ZH5', '8OST', '7XSQ',
#        '4M30', '3U4M', '7DMQ', '8ASW', '6KL9', '6IV9', '2DVI', '7X34',
#        '7ZOL', '5Y7M', '1SJ3', '7KA0', '3WQY', '7NGA', '9IIZ',  # '5H0R', 
#        '4AL6', '3WBM', '2NUE', '7N0C', '1RPU', '4ZLD', '2DU4', '6E9F',
#        '7YND', '6VRB', '9ENC', '2ZH7', '4TUW', '6NMC', '8BGU', '1DFU',
#        '9DTT', '2F8K', '8WM4', '4NH6', '2DU5', '7W0D', '6H66', '2PLY',
#        '7VTN', '8GNA', '5JC7', '7DLZ', '7M5O', '2EZ6', '7YN9', '5HC9',
#        '6WW6', '8HKE', '8EEY', '4KR9', '4X0A', '7Y9X', '1MMS', '2ZKO',
#        '2DER', '1QTQ', '8Z8X', '7QR3', '6SX2', '1VBY', '7Y84', '1ZDK',
#        '8DG7', '5CCB', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ', '8OPT',
#        '2R8S', '2RE8', '1QA6', '6DU4', '8D9I', '8S9T', '5W1H', '6G1X',
#        '2ANN', '5W0O', '3A2K', '2AKE', '4WC6', '4WC2', '4ZT9', '7X7R',
#        '4QOZ', '7RE2', '2B2G', '8OPS', '6XJQ', '3AM1', '2GIC', '6XJW',
#        '4YVI', '1EXD', '8IDF', '1SDS', '2ZH1', '8D1V', '2ZH4', '7YNC',
#        '6AAX', '3KTW', '1P6V', '4RWO', '2BQ5', '8S9U', '1O0B']

# # '6TZ2',

# protein_rna_orphan_ids = ['8Z90', '3OV7', '2UWM', '4J1G', '8D9H', '3OVA', '2ZH9', '6XH1',
#        '4AL5', '7OS0', '7WKP', '2PJP', '4RWP', '4JGN', '4M4O',
#        '7Y80', '4BHH', '6SY6', '8E28', '7P0V', '5WLH', '6DTD', '2ZH2',
#        '3FTF', '7BKQ', '5HSW', '2DR8', '7Y8Y', '6NMD', '4KJI', '7W0C',
#        '8D9G', '7UOB', '7OG0', '7DTE', '8C4U', '2F8T', '7XSP', '5DE5',
#        '8AXF', '2ZHA', '6TQB', '2F8S', '8D9F', '4NGF', '2DR9', '6HTU',
#        '6P7M', '7Y85', '6XH2', '2NUF', '6O5F', '4PMI', '3R9W', '7NQ4',
#        '6IV6', '2OZB', '7U2A', '7XSS', '2ZHB', '8URB', '6SX0', '7X7A',
#        '4TUX', '2ZH6', '7WAH', '7XT4', '8H69', '8SQ9', '7YNB', '6G1S',
#        '8DVS', '6BJV', '9CPO', '6B14', '1ZBH', '8DGA', '7CTT', '8HIO',
#        '7ZOQ', '4FTS', '6BJG', '8H0S', '8JTR', '4IG8', '7Y82', '7NIC',
#        '7UOE', '5DV7', '7X8A', '4KR6', '7R97', '7OZQ', '1JBT', '3PU4',
#        '1VFG', '4JK0', '7KRP', '7ED5', '1S03', '8D4B', '8R3L', '7N0B',
#        '4RWN', '5F5F', '7DOK', '2ZH8', '8WMI', '7SWF', '5MSF', '4QIL',
#        '8WMC', '2DR7', '8D97', '3BOY', '8IUO', '2W2H', '7Y81', '4OOG',
#        '8S9V', '6NMA', '3HTX', '8GU6', '4J39', '5ZW4', '8D8N', '4QI2',
#        '8AF0', '8G6R', '4QQB', '6XH3', '3OUY', '2DRA', '3VJR', '7YNA',
#        '7ORN', '6ZLC', '8SSW', '1EC6', '8E29', '7VG2', '7Y8T', '4L8H',
#        '8T5S', '6GKH', '9IIY', '1FEU', '6CMN', '5F5H', '6W62', '3SIU',
#        '3NVI', '8E40', '8WML', '4M2Z', '3OVS', '7XC7', '8GS2', '9IJ2',
#        '7VG3', '8S9X', '4JNG', '1A9N', '4GHL', '2DR5', '7UO7', '6BJY',
#        '4C7O', '6BJH', '9CMP', '2ZH3', '8EEX', '7XW2', '7XSR', '4PDB',
#        '4S3N', '8D9E', '6H61', '7Y9Y', '7UO4', '2DRB', '2ZH5', '7XSQ',
#        '4M30', '6KL9', '2DVI', '7ZOL', '7NGA', '9IIZ', '4AL6',  # '5H0R',
#        '3WBM', '2NUE', '7N0C', '1RPU', '4ZLD', '7YND', '6VRB', '2ZH7',
#        '4TUW', '6NMC', '1DFU', '2F8K', '8WM4', '4NH6', '6H66', '2PLY',
#        '7VTN', '8GNA', '5JC7', '7M5O', '2EZ6', '7YN9', '6WW6', '8HKE',
#        '8EEY', '4KR9', '7Y9X', '2ZKO', '8Z8X', '6SX2', '7Y84', '1ZDK',
#        '8DG7', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ', '6DU4', '8D9I',
#        '8S9T', '5W1H', '6G1X', '2ANN', '5W0O', '7X7R', '4QOZ', '7RE2',
#        '2B2G', '2GIC', '8IDF', '1SDS', '2ZH1', '8D1V', '2ZH4', '7YNC',
#        '6AAX', '4RWO', '2BQ5', '8S9U']

# protein_rna_non_orphan_ids = ['9ENB', '5VW1', '1GTS', '3NDB', '5NPM', '1EIY', '6E9E', '2HVY',
#        '5UD5', '1VBX', '2RD2', '1SJ4', '7PMQ', '3U56', '7D7V', '3AMU',
#        '8HMY', '3WC2', '1KUQ', '3V7E', '4QG3', '1J1U', '2DET', '8DFV',
#        '1ZJW', '8OMR', '8ISS', '2DU6', '3IVK', '1VC6', '8AW3', '7W0E',
#        '6U8D', '4XCO', '1LNG', '5WT3', '2IY5', '4YVK', '4P3E', '1QRS',
#        '8FZR', '3EPJ', '4PRF', '1EUY', '4WC4', '5DAR', '8YE6', '9ENF',
#        '8PTZ', '7DWH', '5V6X', '3IWN', '7KAB', '1VC0', '4WC3', '1QRU',
#        '1IVS', '6VRC', '1GAX', '7K9M', '9ENE', '1O0C', '3K0J', '3UMY',
#        '6XJZ', '2OJ3', '6DB8', '6B3K', '7MLX', '8CTI', '6DB9', '1VBZ',
#        '2ZZN', '2ZNI', '6ZDQ', '4YYE', '4WC7', '1GTR', '2V3C', '2ZXU',
#        '1QRT', '5D6G', '8PTX', '7QR4', '3AL0', '4YCO', '6V5B', '4WJ4',
#        '4X0B', '1U0B', '4YB1', '2OIH', '3R9X', '4WC5', '6XJY', '1UN6',
#        '1SJF', '4YVJ', '2DR2', '7SLP', '3TUP', '1EUQ', '3FOZ', '8OST',
#        '3U4M', '7DMQ', '8ASW', '6IV9', '7X34', '5Y7M', '1SJ3', '7KA0',
#        '3WQY', '2DU4', '6E9F', '9ENC', '8BGU', '9DTT', '2DU5', '7W0D',
#        '7DLZ', '5HC9', '4X0A', '1MMS', '2DER', '1QTQ', '7QR3', '1VBY',
#        '5CCB', '8OPT', '2R8S', '2RE8', '1QA6', '3A2K', '2AKE', '4WC6',
#        '4WC2', '4ZT9', '8OPS', '6XJQ', '3AM1', '6XJW', '4YVI', '1EXD',
#        '3KTW', '1P6V', '1O0B']

# protein_rna_1v1_ids = ['7XW2', '8FZR', '8DVS', '7Y80', '8D97', '8WMI', '7YN9', '8E28', '7X34', '8E29', '8JTR', '8WM4', '9DTT', '7X7A', '7WKP', '8IDF', '7QR4', '7P0V', '6VRC', '6IV6', '6W62', '6E9E', '4QG3', '1QRT', '2ZH9', '2ZH1', '4J39', '6TQB', '4WC2', '1VC6', '1QRS', '2ZH6', '1O0C', '2DR5', '1LNG', '1VC0', '1O0B', '4ZLD', '2RE8', '6DU4', '1VBZ', '2RD2', '1KUQ', '1QTQ', '2ZH8', '1EXD', '3OVA', '1U0B', '2ZH2', '4YB1', '1VBX', '2DR9', '1SJ4', '4M4O', '2ZH5', '1EUY', '1SJ3', '3R9W', '2DR7', '2ZH4', '4NH6', '5NPM', '2DR8', '1VBY', '2ZH3', '4PDB', '1ZJW', '4PRF', '6ZDQ', '1J1U', '7D7V', '4AL5', '1SJF', '6XH1', '2PJP', '2DRB', '2ZHA', '2OIH', '2ANN', '6CMN', '4AL6', '6XH3', '4JK0', '1GTS', '5ZW4', '2F8K', '4JGN', '2DRA', '6XH2', '5WLH', '1GTR', '5W1H', '3U56', '2ZH7', '1QRU', '5D6G', '4WC6', '1EUQ', '4WC7', '3AMU', '2ZHB', '2DET', '3UMY', '2OJ3', '2DR2', '4WC5', '4WC4', '4WC3', '3TUP', '2AKE', '3U4M', '4X0A', '5WT3', '4X0B']

# protein_rna_multichain_ids = ['9ENB', '5VW1', '8Z90', '3OV7', '2UWM', '3NDB', '4J1G', '8D9H', '7OS0', '1EIY', '2HVY', '4RWP', '5UD5', '4BHH', '6SY6', '3FTF', '7BKQ', '7PMQ', '7Y8Y', '8HMY', '6NMD', '4KJI', '7W0C', '8D9G', '7UOB', '7OG0', '7DTE', '8C4U', '3WC2', '2F8T', '3V7E', '7XSP', '5DE5', '8AXF', '8DFV', '2F8S', '8OMR', '8ISS', '2DU6', '8D9F', '4NGF', '6HTU', '6P7M', '3IVK', '7Y85', '8AW3', '7W0E', '6U8D', '2NUF', '6O5F', '4PMI', '4XCO', '7NQ4', '2OZB', '7U2A', '7XSS', '8URB', '6SX0', '2IY5', '4YVK', '4TUX', '4P3E', '7WAH', '7XT4', '8H69', '8SQ9', '7YNB', '6G1S', '6BJV', '9CPO', '6B14', '1ZBH', '8DGA', '7CTT', '8HIO', '7ZOQ', '4FTS', '6BJG', '8H0S', '3EPJ', '4IG8', '7Y82', '7NIC', '7UOE', '5DV7', '7X8A', '4KR6', '5DAR', '7R97', '7OZQ', '8YE6', '1JBT', '3PU4', '1VFG', '9ENF', '7KRP', '7ED5', '1S03', '8PTZ', '7DWH', '5V6X', '8D4B', '3IWN', '8R3L', '7N0B', '4RWN', '5F5F', '7KAB', '7DOK', '7SWF', '5MSF', '1IVS', '4QIL', '8WMC', '1GAX', '3BOY', '8IUO', '2W2H', '7K9M', '9ENE', '3K0J', '7Y81', '4OOG', '8S9V', '6NMA', '6XJZ', '3HTX', '6DB8', '8GU6', '6B3K', '7MLX', '8CTI', '6DB9', '8D8N', '4QI2', '2ZZN', '8AF0', '8G6R', '2ZNI', '4QQB', '4YYE', '3OUY', '2V3C', '3VJR', '2ZXU', '7YNA', '7ORN', '8PTX', '6ZLC', '8SSW', '1EC6', '7VG2', '7Y8T', '4L8H', '8T5S', '6GKH', '9IIY', '1FEU', '3AL0', '4YCO', '6V5B', '5F5H', '3SIU', '3NVI', '8E40', '8WML', '4M2Z', '3OVS', '7XC7', '8GS2', '9IJ2', '7VG3', '8S9X', '4JNG', '1A9N', '4GHL', '3R9X', '7UO7', '6BJY', '4C7O', '6BJH', '9CMP', '8EEX', '7XSR', '4S3N', '8D9E', '6XJY', '1UN6', '4YVJ', '7SLP', '6H61', '7Y9Y', '7UO4', '3FOZ', '8OST', '7XSQ', '4M30', '7DMQ', '8ASW', '6KL9', '7ZOL', '5Y7M', '7KA0', '3WQY', '7NGA', '9IIZ', '3WBM', '2NUE', '7N0C', '1RPU', '2DU4', '6E9F', '7YND', '6VRB', '9ENC', '4TUW', '6NMC', '8BGU', '1DFU', '2DU5', '7W0D', '6H66', '2PLY', '7VTN', '8GNA', '5JC7', '7DLZ', '2EZ6', '5HC9', '6WW6', '8HKE', '8EEY', '4KR9', '7Y9X', '1MMS', '2ZKO', '2DER', '8Z8X', '7QR3', '6SX2', '7Y84', '1ZDK', '8DG7', '5CCB', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ', '8OPT', '2R8S', '1QA6', '8D9I', '8S9T', '6G1X', '5W0O', '3A2K', '4ZT9', '7X7R', '4QOZ', '7RE2', '2B2G', '8OPS', '6XJQ', '2GIC', '6XJW', '4YVI', '1SDS', '8D1V', '7YNC', '6AAX', '3KTW', '1P6V', '4RWO', '2BQ5', '8S9U']

# # '6TZ2',

# protein_rna_b2021_ids = ['5VW1', '1GTS', '3OV7', '2UWM', '3NDB', '4J1G', '3OVA', '2ZH9',
#        '6XH1', '5NPM', '4AL5', '1EIY', '6E9E', '2PJP', '2HVY', '4RWP',
#        '4JGN', '4M4O', '5UD5', '4BHH', '1VBX', '6SY6', '5WLH',
#        '6DTD', '2ZH2', '3FTF', '2RD2', '1SJ4', '5HSW', '3U56', '7D7V',
#        '2DR8', '3AMU', '6NMD', '4KJI', '3WC2', '2F8T', '1KUQ', '3V7E',
#        '5DE5', '4QG3', '1J1U', '2ZHA', '2DET', '6TQB', '2F8S', '1ZJW',
#        '2DU6', '4NGF', '2DR9', '6HTU', '6P7M', '3IVK', '1VC6', '6XH2',
#        '6U8D', '2NUF', '6O5F', '4PMI', '4XCO', '3R9W', '7NQ4', '1LNG',
#        '6IV6', '2OZB', '2ZHB', '5WT3', '6SX0', '2IY5', '4YVK', '4TUX',
#        '4P3E', '2ZH6', '1QRS', '6G1S', '6BJV', '6B14', '1ZBH', '7CTT',
#        '4FTS', '6BJG', '3EPJ', '4PRF', '4IG8', '1EUY', '5DV7', '4WC4',
#        '4KR6', '5DAR', '1JBT', '3PU4', '1VFG', '4JK0', '7KRP', '1S03',
#        '5V6X', '3IWN', '7N0B', '4RWN', '5F5F', '7KAB', '1VC0', '4WC3',
#        '2ZH8', '1QRU', '5MSF', '1IVS', '4QIL', '2DR7', '6VRC', '1GAX',
#        '3BOY', '2W2H', '7K9M', '1O0C', '3K0J', '4OOG', '6NMA', '3UMY',
#        '3HTX', '2OJ3', '6DB8', '6B3K', '7MLX', '4J39', '5ZW4', '6DB9',
#        '1VBZ', '4QI2', '2ZZN', '2ZNI', '6ZDQ', '4QQB', '6XH3', '4YYE',
#        '4WC7', '3OUY', '1GTR', '2V3C', '2DRA', '3VJR', '2ZXU', '1QRT',
#        '5D6G', '6ZLC', '1EC6', '4L8H', '6GKH', '1FEU', '3AL0', '6CMN',
#        '4YCO', '6V5B', '4WJ4', '5F5H', '6W62', '4X0B', '3SIU', '3NVI',
#        '1U0B', '4M2Z', '3OVS', '4YB1', '4JNG', '2OIH', '1A9N', '4GHL',
#        '2DR5', '3R9X', '6BJY', '4C7O', '6BJH', '4WC5', '2ZH3', '4PDB',
#        '4S3N', '1UN6', '1SJF', '4YVJ', '2DR2', '3TUP', '6H61', '1EUQ',
#        '2DRB', '3FOZ', '2ZH5', '4M30', '3U4M', '7DMQ', '6KL9', '6IV9',
#        '2DVI', '5Y7M', '1SJ3', '7KA0', '3WQY', '4AL6', '3WBM',  # '5H0R', 
#        '2NUE', '7N0C', '1RPU', '4ZLD', '2DU4', '6E9F', '6VRB', '2ZH7',
#        '4TUW', '6NMC', '1DFU', '2F8K', '4NH6', '2DU5', '6H66', '2PLY',
#        '5JC7', '7M5O', '2EZ6', '5HC9', '4KR9', '4X0A', '1MMS', '2ZKO',
#        '2DER', '1QTQ', '6SX2', '1VBY', '1ZDK', '5CCB', '4X4O', '6CF2',
#        '6KLB', '1OOA', '5UDZ', '2R8S', '2RE8', '1QA6', '6DU4', '5W1H',
#        '6G1X', '2ANN', '5W0O', '3A2K', '2AKE', '4WC6', '4WC2', '4ZT9',
#        '4QOZ', '2B2G', '3AM1', '2GIC', '4YVI', '1EXD', '1SDS', '2ZH1',
#        '2ZH4', '6AAX', '3KTW', '1P6V', '4RWO', '2BQ5', '1O0B']

# # '6TZ2',

# protein_rna_b2021_orphan_ids = ['3OV7', '2UWM', '4J1G', '3OVA', '2ZH9', '6XH1', '4AL5', '2PJP',
#        '4RWP', '4JGN', '4M4O', '4BHH', '6SY6', '5WLH', '6DTD',
#        '2ZH2', '3FTF', '5HSW', '2DR8', '6NMD', '4KJI', '2F8T', '5DE5',
#        '2ZHA', '6TQB', '2F8S', '4NGF', '2DR9', '6HTU', '6P7M', '6XH2',
#        '2NUF', '6O5F', '4PMI', '3R9W', '7NQ4', '6IV6', '2OZB', '2ZHB',
#        '6SX0', '4TUX', '2ZH6', '6G1S', '6BJV', '6B14', '1ZBH', '7CTT',
#        '4FTS', '6BJG', '4IG8', '5DV7', '4KR6', '1JBT', '3PU4', '1VFG',
#        '4JK0', '7KRP', '1S03', '7N0B', '4RWN', '5F5F', '2ZH8', '5MSF',
#        '4QIL', '2DR7', '3BOY', '2W2H', '4OOG', '6NMA', '3HTX', '4J39',
#        '5ZW4', '4QI2', '4QQB', '6XH3', '3OUY', '2DRA', '3VJR', '6ZLC',
#        '1EC6', '4L8H', '6GKH', '1FEU', '6CMN', '5F5H', '6W62', '3SIU',
#        '3NVI', '4M2Z', '3OVS', '4JNG', '1A9N', '4GHL', '2DR5', '6BJY',
#        '4C7O', '6BJH', '2ZH3', '4PDB', '4S3N', '6H61', '2DRB', '2ZH5',
#        '4M30', '6KL9', '2DVI', '4AL6', '3WBM', '2NUE', '7N0C',  # '5H0R',
#        '1RPU', '4ZLD', '6VRB', '2ZH7', '4TUW', '6NMC', '1DFU', '2F8K',
#        '4NH6', '6H66', '2PLY', '5JC7', '7M5O', '2EZ6', '4KR9', '2ZKO',
#        '6SX2', '1ZDK', '4X4O', '6CF2', '6KLB', '1OOA', '5UDZ', '6DU4',
#        '5W1H', '6G1X', '2ANN', '5W0O', '4QOZ', '2B2G', '2GIC', '1SDS',
#        '2ZH1', '2ZH4', '6AAX', '4RWO', '2BQ5']

# protein_rna_b2021_non_orphan_ids = ['5VW1', '1GTS', '3NDB', '5NPM', '1EIY', '6E9E', '2HVY', '5UD5',
#        '1VBX', '2RD2', '1SJ4', '3U56', '7D7V', '3AMU', '3WC2', '1KUQ',
#        '3V7E', '4QG3', '1J1U', '2DET', '1ZJW', '2DU6', '3IVK', '1VC6',
#        '6U8D', '4XCO', '1LNG', '5WT3', '2IY5', '4YVK', '4P3E', '1QRS',
#        '3EPJ', '4PRF', '1EUY', '4WC4', '5DAR', '5V6X', '3IWN', '7KAB',
#        '1VC0', '4WC3', '1QRU', '1IVS', '6VRC', '1GAX', '7K9M', '1O0C',
#        '3K0J', '3UMY', '2OJ3', '6DB8', '6B3K', '7MLX', '6DB9', '1VBZ',
#        '2ZZN', '2ZNI', '6ZDQ', '4YYE', '4WC7', '1GTR', '2V3C', '2ZXU',
#        '1QRT', '5D6G', '3AL0', '4YCO', '6V5B', '4WJ4', '4X0B', '1U0B',
#        '4YB1', '2OIH', '3R9X', '4WC5', '1UN6', '1SJF', '4YVJ', '2DR2',
#        '3TUP', '1EUQ', '3FOZ', '3U4M', '7DMQ', '6IV9', '5Y7M', '1SJ3',
#        '7KA0', '3WQY', '2DU4', '6E9F', '2DU5', '5HC9', '4X0A', '1MMS',
#        '2DER', '1QTQ', '1VBY', '5CCB', '2R8S', '2RE8', '1QA6', '3A2K',
#        '2AKE', '4WC6', '4WC2', '4ZT9', '3AM1', '4YVI', '1EXD', '3KTW',
#        '1P6V', '1O0B']

# protein_rna_a2021_ids = ['9ENB', '8Z90', '8D9H', '7OS0', '7WKP', '7Y80', '8E28', '7P0V',
#        '7BKQ', '7PMQ', '7Y8Y', '8HMY', '7W0C', '8D9G', '7UOB', '7OG0',
#        '7DTE', '8C4U', '7XSP', '8AXF', '8DFV', '8OMR', '8ISS', '8D9F',
#        '7Y85', '8AW3', '7W0E', '7U2A', '7XSS', '8URB', '7X7A', '7WAH',
#        '7XT4', '8H69', '8SQ9', '7YNB', '8DVS', '9CPO', '8FZR', '8DGA',
#        '8HIO', '7ZOQ', '8H0S', '8JTR', '7Y82', '7NIC', '7UOE', '7X8A',
#        '7R97', '7OZQ', '8YE6', '9ENF', '7ED5', '8PTZ', '7DWH', '8D4B',
#        '8R3L', '7DOK', '8WMI', '7SWF', '8WMC', '8D97', '8IUO', '9ENE',
#        '7Y81', '8S9V', '6XJZ', '8GU6', '8CTI', '8D8N', '8AF0', '8G6R',
#        '7YNA', '7ORN', '8PTX', '8SSW', '8E29', '7VG2', '7Y8T', '8T5S',
#        '7QR4', '9IIY', '8E40', '8WML', '7XC7', '8GS2', '9IJ2', '7VG3',
#        '8S9X', '7UO7', '9CMP', '8EEX', '7XW2', '7XSR', '8D9E', '6XJY',
#        '7SLP', '7Y9Y', '7UO4', '8OST', '7XSQ', '8ASW', '7X34', '7ZOL',
#        '7NGA', '9IIZ', '7YND', '9ENC', '8BGU', '9DTT', '8WM4', '7W0D',
#        '7VTN', '8GNA', '7DLZ', '7YN9', '6WW6', '8HKE', '8EEY', '7Y9X',
#        '8Z8X', '7QR3', '7Y84', '8DG7', '8OPT', '8D9I', '8S9T', '7X7R',
#        '7RE2', '8OPS', '6XJQ', '6XJW', '8IDF', '8D1V', '7YNC', '8S9U']

# protein_rna_a2021_orphan_ids = ['8Z90', '8D9H', '7OS0', '7WKP', '7Y80', '8E28', '7P0V', '7BKQ',
#        '7Y8Y', '7W0C', '8D9G', '7UOB', '7OG0', '7DTE', '8C4U', '7XSP',
#        '8AXF', '8D9F', '7Y85', '7U2A', '7XSS', '8URB', '7X7A', '7WAH',
#        '7XT4', '8H69', '8SQ9', '7YNB', '8DVS', '9CPO', '8DGA', '8HIO',
#        '7ZOQ', '8H0S', '8JTR', '7Y82', '7NIC', '7UOE', '7X8A', '7R97',
#        '7OZQ', '7ED5', '8D4B', '8R3L', '7DOK', '8WMI', '7SWF', '8WMC',
#        '8D97', '8IUO', '7Y81', '8S9V', '8GU6', '8D8N', '8AF0', '8G6R',
#        '7YNA', '7ORN', '8SSW', '8E29', '7VG2', '7Y8T', '8T5S', '9IIY',
#        '8E40', '8WML', '7XC7', '8GS2', '9IJ2', '7VG3', '8S9X', '7UO7',
#        '9CMP', '8EEX', '7XW2', '7XSR', '8D9E', '7Y9Y', '7UO4', '7XSQ',
#        '7ZOL', '7NGA', '9IIZ', '7YND', '8WM4', '7VTN', '8GNA', '7YN9',
#        '6WW6', '8HKE', '8EEY', '7Y9X', '8Z8X', '7Y84', '8DG7', '8D9I',
#        '8S9T', '7X7R', '7RE2', '8IDF', '8D1V', '7YNC', '8S9U']

# protein_rna_a2021_non_orphan_ids = ['9ENB', '7PMQ', '8HMY', '8DFV', '8OMR', '8ISS', '8AW3', '7W0E',
#        '8FZR', '8YE6', '9ENF', '8PTZ', '7DWH', '9ENE', '6XJZ', '8CTI',
#        '8PTX', '7QR4', '6XJY', '7SLP', '8OST', '8ASW', '7X34', '9ENC',
#        '8BGU', '9DTT', '7W0D', '7DLZ', '7QR3', '8OPT', '8OPS', '6XJQ',
#        '6XJW']

val_ids = ['3BWP', '255D', '6E80', '4FAQ', '3Q50', '2ZY6', '4E8V', '7KD1', '3AM1', '2Q1R', '3GCA', '3ND3', '3DHS', '3NPN', '7M5O', '4E8P', '4E8Q', '4RBQ', '1U9S', '4WJ4', '6UES', '4EN5', '4E8M', '4C40', '6TF3', '5C5W', '4CS1', '4E8N', '5DA6', '6TB7', '4P8Z', '2A2E', '6IV9', '2A64', '5HSW', '413D', '3R4F', '2DVI', '4GMA', '6TFE', '3D0M', '4DS6', '387D', '7D7W', '6TF1', '6UET', '6T3S', '6DTD', '6PQ7', '4AOB']

ids_prefix_suffix_mapping = {
    'all': [
        (all_data_ids, 'All', '.csv'),
        (all_data_orphan_ids, 'All', '_orphan.csv'),
        (all_data_non_orphan_ids, 'All', '_non_orphan.csv'),
        (all_data_b2021_ids, 'All_b2021', '.csv'),
        (all_data_b2021_orphan_ids, 'All_b2021', '_orphan.csv'),
        (all_data_b2021_non_orphan_ids, 'All_b2021', '_non_orphan.csv'),
        (all_data_a2021_ids, 'All_a2021', '.csv'),
        (all_data_a2021_orphan_ids, 'All_a2021', '_orphan.csv'),
        (all_data_a2021_non_orphan_ids, 'All_a2021', '_non_orphan.csv'),
    ],
    'protein_rna': [
        (protein_rna_ids, 'RNA-Protein', '.csv'),
        (protein_rna_orphan_ids, 'RNA-Protein', '_orphan.csv'),
        (protein_rna_non_orphan_ids, 'RNA-Protein', '_non_orphan.csv'),
        (protein_rna_b2021_ids, 'RNA-Protein_b2021', '.csv'),
        (protein_rna_b2021_orphan_ids, 'RNA-Protein_b2021', '_orphan.csv'),
        (protein_rna_b2021_non_orphan_ids, 'RNA-Protein_b2021', '_non_orphan.csv'),
        (protein_rna_a2021_ids, 'RNA-Protein_a2021', '.csv'),
        (protein_rna_a2021_orphan_ids, 'RNA-Protein_a2021', '_orphan.csv'),
        (protein_rna_a2021_non_orphan_ids, 'RNA-Protein_a2021', '_non_orphan.csv'),
        (protein_rna_1v1_ids, 'RNA-Protein', '_1v1.csv'),
        (protein_rna_multichain_ids, 'RNA-Protein', '_multichain.csv'),
    ],
    'rna': [
        (rna_monomers_ids, 'RNA_Monomers', '.csv'),
        (rna_monomers_orphan_ids, 'RNA_Monomers', '_orphan.csv'),
        (rna_monomers_non_orphan_ids, 'RNA_Monomers', '_non_orphan.csv'),
        (rna_monomers_b2021_ids, 'RNA_Monomers_b2021', '.csv'),
        (rna_monomers_b2021_orphan_ids, 'RNA_Monomers_b2021', '_orphan.csv'),
        (rna_monomers_b2021_non_orphan_ids, 'RNA_Monomers_b2021', '_non_orphan.csv'),
        (rna_monomers_a2021_ids, 'RNA_Monomers_a2021', '.csv'),
        (rna_monomers_a2021_orphan_ids, 'RNA_Monomers_a2021', '_orphan.csv'),
        (rna_monomers_a2021_non_orphan_ids, 'RNA_Monomers_a2021', '_non_orphan.csv'),
    ],
}

# full_columns = ['Id', 'exp_db_id', 'FileName', 'ProteinSequence', 'ProteinChainIDs',
#        'ProteinLength', 'NumberProteins', 'AAMotif', 'AAC', 'AAPproteinRNA',
#        'RNASequence', 'RNAChainIDs', 'RNALength', 'NumberRNAs', 'RNAMotif',
#        'RNAMotifLength', 'RNAModified', 'ContactList',
#        'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
#        'ProteinRNAInterfaceArea', 'ProteinRNAInterfaceRatio', 'Free_energy',
#        'Binding_affinity_kd', 'Seed', 'Model', 'af3_protein_pTM',
#        'af3_rna_pTM', 'af3_protein_ipTM', 'af3_rna_ipTM',
#        'af3_protein_pLDDT_avg', 'af3_rna_pLDDT_avg', 'af3_global_pae_avg',
#        'af3_chain_pair_pae_min', 'af3_fraction_disordered', 'af3_has_clash',
#        'af3_ranking_score', 'af3_protein_MSA', 'af3_rna_MSA',
#        'rna_prot_interface_atom_ids', 'interface_rna_atom_ids', 'Complex_RMSD',
#        'Protein_RMSD', 'RNA_RMSD', 'Protein_LDDT', 'RNA_LDDT', 'Protein_TM',
#        'RNA_TM', 'Complex_TM', 'Complex_LDDT', 'RNA_INF', 'RNA_DI',
#        'RNA_f1_score', 'RNA_precision', 'RNA_recall', 'RNA_mcc', 'RNA_wl',
#        'RNA_GlycosidicBond', 'RNA_SugarPucker', 'RNA_GammaAngle',
#        'RNAmotif_score', 'RNA_Stems', 'RNA_HairpinLoops', 'RNA_InternalLoops',
#        'RNA_MultibranchLoops', 'RNA_DanglingEnds', 'RNA_Pseudoknots',
#        'RNA_isDuplex', 'RNAMotif_isDuplex', 'RNA_LDDT_allAtom',
#        'Protein_LDDT_allAtom', 'Complex_LDDT_allAtom', 'iLDDT', 'QS_global',
#        'Dockq', 'ICS', 'ChainIDpairList_rnaRNA', 'Hbond_rnaRNA',
#        'vdWbond_rnaRNA', 'RNArnaInterfaceArea', 'RNArnaInterfaceRatio',
#        'rna_rna_interface_atom_ids']

full_columns = ['Id', 'exp_db_id', 'FileName', 'ProteinSequence', 'ProteinChainIDs',
       'ProteinLength', 'NumberProteins', 'AAMotif', 'AAC', 'AAPproteinRNA',
       'RNASequence', 'RNAChainIDs', 'RNALength', 'NumberRNAs', 'RNAMotif',
       'RNAMotifLength', 'RNAModified', 'ContactList',
       'ChainIDpairList_proteinRNA', 'Hbond_proteinRNA', 'vdWbond_proteinRNA',
       'ProteinRNAInterfaceArea', 'ProteinRNAInterfaceRatio', 'Free_energy',
       'Binding_affinity_kd', 'Seed', 'Model', 'af3_protein_pTM',
       'af3_rna_pTM', 'af3_protein_ipTM', 'af3_rna_ipTM',
       'af3_protein_pLDDT_avg', 'af3_rna_pLDDT_avg', 'af3_global_pae_avg',
       'af3_chain_pair_pae_min', 'af3_fraction_disordered', 'af3_has_clash',
       'af3_ranking_score', 'af3_protein_MSA', 'af3_rna_MSA',
       'rna_prot_interface_atom_ids', 'interface_rna_atom_ids', 'Complex_RMSD',
       'Protein_RMSD', 'RNA_RMSD', 'Protein_LDDT', 'RNA_LDDT', 'Protein_TM',
       'RNA_TM', 'Complex_TM', 'Complex_LDDT', 'RNA_INF', 'RNA_DI',
       'RNA_f1_score', 'RNA_precision', 'RNA_recall', 'RNA_mcc', 'RNA_wl',
       'RNA_GlycosidicBond', 'RNA_SugarPucker', 'RNA_GammaAngle',
       'RNAmotif_score', 'RNA_Stems', 'RNA_HairpinLoops', 'RNA_InternalLoops',
       'RNA_MultibranchLoops', 'RNA_DanglingEnds', 'RNA_Pseudoknots',
       'RNA_isDuplex', 'RNAMotif_isDuplex', 'RNAMotif_inExp',
       'RNA_ElectrostaticPotential', 'RNA_GC_Content',
       'RNA_SequenceComplexity', 'RNA_Shannon_Entropy_K1',
       'RNA_Shannon_Entropy_K2', 'RNA_GC_Skew', 'RNA_AU_Skew',
       'RNA_AU_Content', 'RNA_MFE_Value', 'RNA_MFE_Structure',
       'RNA_Ensemble_Diversity', 'RNA_Centroid_Energy', 'RNA_msa_size',
       'Protein_msa_size', 'AAmatch_score', 'RNAmatch_score',
       'RNA_LDDT_allAtom', 'Protein_LDDT_allAtom', 'Complex_LDDT_allAtom',
       'iLDDT', 'QS_global', 'Dockq', 'ICS', 'RNAComplex_TM']

def get_csvs(directory, pattern=''):
    return list(Path(directory).glob(f"{pattern}.csv"))

def load_csv_to_df(csv_path, sep=','):
    try:
        return pd.read_csv(csv_path)
    except Exception as e:
        print(e)
        return pd.read_csv(csv_path, sep=';')
    #     try:
    #         return pd.read_csv(csv_path, sep=',')
    #     except Exception as e:
    #         print(e)
    #         return pd.read_csv(csv_path, sep=';')
    # finally:
    #     raise UserWarning('Cannot read csv file for %s', csv_path)

def get_rmsd(x):
    if isinstance(x, str):
        value = ast.literal_eval(x)
        if isinstance(value, list):
            return float(value[0])
        else: 
            return float(value)
    elif isinstance(x, (int, float)):
        if pd.isna(x):
            return np.nan
        else:
            return float(x)
    else:
        return np.nan

def get_crmsd_allatom(x):
    if isinstance(x, str):
        value = ast.literal_eval(x)
        if isinstance(value, list):
            return float(value[2])
        else: 
            return float(value)
    elif isinstance(x, (int, float)):
        if pd.isna(x):
            return np.nan
        else:
            return float(x)
    else:
        return np.nan

def get_max_rna_length(row):
    """
    Extract the maximum RNA length from the 'RNA_msa_size' column.
    The column is expected to be a string representation of a dictionary.
    """
    try:
        rnalengths = ast.literal_eval(str(row['RNALength']))
    except (ValueError, SyntaxError):
        print(f"Error parsing RNALength for row: {row}")
        print(row['RNALength'])
        rnalengths = row['RNALength']
    if isinstance(rnalengths, int):
        return rnalengths
    elif isinstance(rnalengths, str):
        rnalengths = [int(l) for l in rnalengths.split(',')]
    elif isinstance(rnalengths, list):
        rnalengths = [int(l) for l in rnalengths]
    return np.max([int(x) for x in rnalengths])


def main(directory='results/new_eval_iris', algorithm='rnafold'):
    
    initial_data = get_csvs(directory, pattern=f"{algorithm}_rna")
    initial_data += get_csvs(directory, pattern=f"{algorithm}_protein_rna")

    print(initial_data)
    # print([len(pd.read_csv(p)) for p in initial_data])
    
    df = pd.concat([load_csv_to_df(csv) for csv in initial_data]).reindex(columns=full_columns)

    # drop rows where these are NaN
    # df = df.dropna(subset=['RNA_RMSD', 'RNA_LDDT', 'RNA_TM'])

    df['RNA_RMSD'] = df['RNA_RMSD'].apply(get_rmsd)
    df['Complex_RMSD'] = df['Complex_RMSD'].apply(get_rmsd)
    df['Complex_RMSD'] = df.apply(
        lambda x: float(x['Complex_RMSD']) if pd.notna(x['Complex_RMSD']) else float(x['RNA_RMSD']),
        axis=1
    )
    df['Complex_RMSD_allAtom'] = df['Complex_RMSD'].apply(get_crmsd_allatom)
    df['iLDDT'] = df['iLDDT'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['RNA_LDDT'] = df['RNA_LDDT'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['RNA_TM'] = df['RNA_TM'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['Complex_LDDT'] = df['Complex_LDDT'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['Complex_LDDT'] = df.apply(
        lambda x: float(x['Complex_LDDT']) if pd.notna(x['Complex_LDDT']) else float(x['RNA_LDDT']),
        axis=1
    )
    df['Complex_TM'] = df['Complex_TM'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['Complex_TM'] = df.apply(
        lambda x: float(x['Complex_TM']) if pd.notna(x['Complex_TM']) else float(x['RNA_TM']),
        axis=1
    )

    df['QS_global'] = df['QS_global'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['Dockq'] = df['Dockq'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['ICS'] = df['ICS'].apply(lambda x: float(x) if pd.notna(x) else np.nan)

    df['RNAComplex_TM'] = df['RNAComplex_TM'].apply(lambda x: float(x) if pd.notna(x) else np.nan)
    df['Complex_LDDT_allAtom'] = df['Complex_LDDT_allAtom'].apply(
        lambda x: float(x) if pd.notna(x) else np.nan
    )
    df['RNA_LDDT_allAtom'] = df['RNA_LDDT_allAtom'].apply(
        lambda x: float(x) if pd.notna(x) else np.nan
    )
    df['algorithm'] = algorithm

    df['NumberRNAs'] = df['NumberRNAs'].apply(lambda x: int(x) if pd.notna(x) else 0)
    df['NumberProteins'] = df['NumberProteins'].apply(lambda x: int(x) if pd.notna(x) else 0)
    df['RNALength'] = df.apply(get_max_rna_length, axis=1)

    return df




if __name__ == '__main__':
    import argparse

    # directory = 'results/new_eval_iris_latest'
    directory = 'results/iris_eval_neurips_2026_04_05_2026'
    outdir = 'results/csvs_neurips'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    BASELINE_MODEL = "alphafold"

    model_configs = {
        # 'rnafold': {'required_for_common': True},
        # 'rnaformer': {'required_for_common': True},
        # 'spotrna': {'required_for_common': True},
        
        # 'rnaformerN100': {'required_for_common': True},  # 479 full benchmark
        'rnaformerN100': {'required_for_common': True},  # 479 full benchmark
        'rnaformerN10000': {'required_for_common': True},  # 91 rnaonly
        'rnaformerN100SHSandMSA': {'required_for_common': True},  # 479 full benchmark
        'dssrN100SHSandMSA': {'required_for_common': True},  # 91 rnaonly
        'rnaformerN100ShuffleandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100Shuffle': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100Random': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100RandomandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN5000': {'required_for_common': True},  # 91 rnaonly
        'rnaformerN20000': {'required_for_common': True},  # 91 rnaonly
        'rhofold': {'required_for_common': False},  # 91 rnaonly
        'dssrN100Rhofold': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100DenoiseRhofold': {'required_for_common': False},  # 91 rnaonly
        'spotrnaN100DenoiseRhofold': {'required_for_common': False},  # 91 rnaonly
        'rnafoldN100DenoiseRhofold': {'required_for_common': False},  # 91 rnaonly
        'rnafoldN100Rhofold': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100Rhofold': {'required_for_common': False},  # 91 rnaonly
        'spotrnaN100Rhofold': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100S0': {'required_for_common': False},  # 91 rnaonly
        'dssrN100': {'required_for_common': False},  # 479 benchmark ; 91 rnaonly
        'rnaformerN100S42': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100S137': {'required_for_common': False},  # 91 rnaonly
        'rnafoldN100Denoise': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100Denoise': {'required_for_common': False},  # 91 rnaonly
        'spotrnaN100Denoise': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN20SHSandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN50SHSandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN200SHSandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN500SHSandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN1000SHSandMSA': {'required_for_common': False},  # 91 rnaonly
        'rnaformerN100SHSandMSARhofold': {'required_for_common': False},  # 91 rnaonly
        # 'rnaformerN10000': {'required_for_common': False},  # 479 full benchmark
        # 'rnaformerN5000': {'required_for_common': False},  # 479 full benchmark
        # 'rnaformerN20000': {'required_for_common': False},  # 479 full benchmark
        'alphafold': {'required_for_common': True},  # 479 full benchmark
        'rnafoldN100': {'required_for_common': True},  # 479 full benchmark
        'spotrnaN100': {'required_for_common': True},  # 479 full benchmark
        'noMSA': {'required_for_common': True},  # 479 full benchmark
        # 'rhofoldRNAformerN100': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldRNAfoldN100': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldSpotrnaN100': {'required_for_common': False},  # 91/90 rnaonly
        # 'rhofoldDSSRN100': {'required_for_common': False},  # 91 rnaonly
        # 'rhofold': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldRNAformerN100Denoise': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldRNAfoldN100Denoise': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldSpotrnaN100Denoise': {'required_for_common': False},  # 91 rnaonly
        # 'rhofoldRNAformerN100PlusMSA': {'required_for_common': False},  # 91 rnaonly
        # 'dssrN100': {'required_for_common': False},  # 477 full benchmark or so
        # 'rnaformerN100Denoise': {'required_for_common': False},  # 91 rnaonly
        # 'rnafoldN100Denoise': {'required_for_common': False},  # 91 rnaonly
        # 'spotrnaN100Denoise': {'required_for_common': False},  # 91 rnaonly
    
        # Add new partial runs here
        # 'rnaformerN100_msa_plus': {'required_for_common': False},
        # 'gtN100_msa_plus': {'required_for_common': False},
    }
    
    model_list = list(model_configs.keys())
    
    dfs_by_model = {}
    for model in model_list:
        df_model = main(directory, model)
        dfs_by_model[model] = df_model
        print(model, "loaded rows:", len(df_model), "unique ids:", df_model["exp_db_id"].nunique())
    
    for group_name, mappings in ids_prefix_suffix_mapping.items():
        for ids, prefix, suffix in mappings:
            ids_set = set(ids)
    
            subsets = {}
            required_id_sets = []
    
            for model, df in dfs_by_model.items():
                d = df[df['exp_db_id'].isin(ids_set)].copy()
                subsets[model] = d
    
                if model_configs[model].get('required_for_common', True):
                    if d.empty:
                        print(f"WARNING: required model {model} empty for {prefix}{suffix}")
                    else:
                        required_id_sets.append(set(d['exp_db_id']))
    
            common_ids = set.intersection(*required_id_sets) if required_id_sets else set()

            print(prefix, suffix, "core common ids:", len(common_ids))
    
            if BASELINE_MODEL in subsets:
                baseline_ids = set(subsets[BASELINE_MODEL]['exp_db_id'])
            else:
                baseline_ids = common_ids
    
            for model, d in subsets.items():
                if model_configs[model].get('required_for_common', True):
                    selected_ids = common_ids
                    comparison_mode = "core_common"
                else:
                    # Partial models do not shrink the core common set.
                    # They are written for IDs shared with baseline within this subset.
                    selected_ids = set(d['exp_db_id']) & baseline_ids
                    comparison_mode = "paired_to_baseline"
    
                d_common = (
                    d[d['exp_db_id'].isin(selected_ids)]
                    .reindex(columns=full_columns + ['algorithm'])
                )
    
                dupes = (
                    d_common
                    .groupby('exp_db_id')
                    .size()
                    .loc[lambda s: s > 1]
                )
    
                if not dupes.empty:
                    print("DUPLICATES for", prefix, suffix, model, ":", dupes.index.tolist())
    
                print(
                    prefix, suffix, model,
                    "mode =", comparison_mode,
                    "n_rows =", len(d_common),
                    "n_ids =", d_common["exp_db_id"].nunique(),
                )
    
                out_path = Path(outdir, f"{prefix}_{model}{suffix}")
                d_common.to_csv(out_path, index=False)
    
    # model_list = [
    #     'rnafold',
    #     'rnaformer',
    #     'spotrna',
    #     'rnaformerN100',
    #     'rnaformerN10000',
    #     'rnaformerN5000',
    #     'rnaformerN20000',
    #     'alphafold',
    #     'rnafoldN100',
    #     'spotrnaN100',
    #     'noMSA',
    # ]

    # # 1) load a df for each algorithm
    # dfs_by_model = {}
    # for model in model_list:
    #     dfs_by_model[model] = main(directory, model)

    # # 2) for each (prefix, suffix) compute common IDs and write per model
    # for group_name, mappings in ids_prefix_suffix_mapping.items():
    #     for ids, prefix, suffix in mappings:
    #         ids_set = set(ids)

    #         # per-model subsets and the sets of exp_db_id actually present
    #         subsets = {}
    #         id_sets = []

    #         for model, df in dfs_by_model.items():
    #             d = df[df['exp_db_id'].isin(ids_set)]
    #             subsets[model] = d
    #             id_sets.append(set(d['exp_db_id']))

    #         # minimal common set of IDs across all algorithms
    #         if id_sets:
    #             common_ids = set.intersection(*id_sets)
    #         else:
    #             common_ids = set()

    #         # write one CSV per algorithm, restricted to common_ids
    #         for model, d in subsets.items():
    #             d_common = d[d['exp_db_id'].isin(common_ids)].reindex(columns=full_columns + ['algorithm'])
    #             # right after you build d_common
    #             dupes = (
    #                 d_common
    #                 .groupby('exp_db_id')
    #                 .size()
    #                 .loc[lambda s: s > 1]
    #             )
                
    #             if not dupes.empty:
    #                 print("DUPLICATES for", prefix, suffix, model, ":", dupes.index.tolist())

    #             print(prefix, suffix, model, 'n_ids =', len(d_common))
    #             out_path = Path(outdir, f"{prefix}_{model}{suffix}")
    #             d_common.to_csv(out_path, index=False)
