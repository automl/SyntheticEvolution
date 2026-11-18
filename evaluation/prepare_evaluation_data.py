import shutil
import os
import urllib.request
import argparse

from pathlib import Path

# All Ids of RNA and RNA-Protein samples from PDB (not all predicted ones, just all at reasonable resolution.)
a2021_rp = ['6fql', '6nbj', '6ww6', '6wxq', '6xjq', '6xjw', '6xjy', '6xjz', '6xn3', '6xn4', '6xn5', '6xn7', '7bah', '7bge', '7bkp', '7bkq', '7d8o', '7dfg', '7dfh', '7dlz', '7doi', '7dok', '7dte', '7dwh', '7ed5', '7eiu', '7exa', '7f1m', '7f3i', '7f3j', '7f3k', '7f3l', '7lj3', '7nga', '7nha', '7nhc', '7nhx', '7ni0', '7nic', '7niq', '7nj3', '7nk1', '7nka', '7nkc', '7oea', '7oeb', '7og0', '7ojk', '7ojl', '7ojn', '7orm', '7orn', '7oro', '7os0', '7ozq', '7ozr', '7p0v', '7pdv', '7pmq', '7qdy', '7qdz', '7qr3',
'7qr4', '7qtl', '7r97', '7r9g', '7rdx', '7rdz', '7re0', '7re2', '7rgu', '7rzz', '7s02', '7slp', '7sop', '7sos', '7sov', '7sow', '7sva',
'7swf', '7tj2', '7to0', '7tqv', '7tuv', '7u2a', '7uj1', '7umk', '7uml', '7uo0', '7uo1', '7uo4', '7uo5', '7uo7', '7uob', '7uoe',
'7uu3', '7uu4', '7uu5', '7uzw', '7uzx', '7uzz', '7v00', '7v2z', '7v4g', '7v6c', '7vg2', '7vg3', '7vkl', '7vsj', '7vtn', '7w0c',
'7w0d', '7w0e', '7w9s', '7wah', '7wkp', '7wkv', '7wl0', '7wnu', '7x34', '7x7a', '7x7r', '7x8a', '7xc7', '7xd8', '7xpl',
'7xs4', '7xsp', '7xsq', '7xsr', '7xss', '7xt4', '7xw2', '7xwz', '7y7p', '7y7q', '7y80', '7y81', '7y82', '7y84', '7y85', '7y8t',
'7y8y', '7y9x', '7y9y', '7yew', '7yex', '7yey', '7yfq', '7yfy', '7ygh', '7yn9', '7yna', '7ynb', '7ync', '7ynd', '7ypw', '7yr7',
'7yr8', '7z43', '7z4o', '7z52', '7zgp', '7zgr', '7zhh', '7zlq', '7zol', '7zoq', '8acb', '8acc', '8af0', '8ane', '8as6', '8asw',
'8aw3', '8axb', '8axf', '8b0r', '8b6z', '8b9g', '8b9i', '8bao', '8bgu', '8bh7', '8bh8', '8bmw', '8bvm', '8c4t', '8c4u',
'8c4v', '8cbw', '8ci5', '8cti', '8cx2', '8d1v', '8d29', '8d4b', '8d8n', '8d97', '8d9e', '8d9f', '8d9g', '8d9h', '8d9i', '8dfo',
'8dfv', '8dg7', '8dga', '8dk7', '8do6', '8dp3', '8dvs', '8e28', '8e29', '8e3i', '8e40', '8e5l', '8e5p', '8e6w', '8edj', '8eex',
'8eey', '8enk', '8ey6', '8ey7', '8ey8', '8fd2', '8fn6', '8fnc', '8fnf', '8ftm', '8fuk', '8fvi', '8fzr', '8g6r', '8g9s', '8g9t', '8gaf',
'8gna', '8gs2', '8gu6', '8gw1', '8gwe', '8gwf', '8gwg', '8gwi', '8gwk', '8gwm', '8gwn', '8gwo', '8h0i', '8h0s', '8h1b',
'8h69', '8h8e', '8h9d', '8hio', '8hke', '8hmy', '8id2', '8idf', '8ipm', '8iss', '8iuo', '8j62', '8j72', '8jh1', '8jsn', '8jtr', '8jy0',
'8k0k', '8k34', '8kag', '8kd9', '8kda', '8ki8', '8omr', '8ops', '8opt', '8ost', '8p0b', '8p0g', '8p0u', '8p1m', '8p1n', '8pdl',
'8pdm', '8pdo', '8pdp', '8pdq', '8pdr', '8pds', '8pnf', '8poh', '8psq', '8psu', '8pt7', '8ptj', '8ptx', '8ptz', '8pu0', '8q3z',
'8q9t', '8qca', '8qgt', '8qjk', '8r3l', '8r3z', '8r65', '8r6u', '8r8r', '8ran', '8rao', '8rn1', '8s9t', '8s9u', '8s9v', '8s9x', '8sj7',
'8snx', '8sny', '8sq9', '8ssw', '8sxu', '8t29', '8t2a', '8t2b', '8t2o', '8t5s', '8t65', '8t66', '8tdv', '8tdw', '8thq', '8ube',
'8ud3', '8ud5', '8urb', '8viv', '8vm8', '8vm9', '8vma', '8vu0', '8wce', '8wcs', '8wm4', '8wmc', '8wmi', '8wml', '8wzc',
'8x0n', '8x5d', '8xko', '8xpo', '8ye6', '8yh9', '8yhd', '8yhe', '8z4l', '8z85', '8z8j', '8z8x', '8z90', '8z97', '8z98', '9asi',
'9aur', '9aus', '9avr', '9c68', '9c69', '9cmp', '9cpo', '9dn4', '9dtt', '9enb', '9enc', '9ene', '9enf', '9eyj', '9f2r', '9f37',
'9fvd', '9iiy', '9iiz', '9ij2', '9ij3', '9ij5']


b2021_rp = ['1a9n', '1aq3', '1aq4', '1av6', '1b7f', '1bmv', '1c9s', '1cvj', '1dfu', '1di2', '1ec6', '1eiy', '1euq',
'1euy', '1exd', '1f8v', '1feu', '1fjg', '1fxl', '1g1x', '1g2e', '1g59', '1gax', '1gtn', '1gtr', '1gts', '1hnx', '1hnz', '1hr0', '1i94',
'1i95', '1i96', '1i97', '1ibk', '1ibl', '1ibm', '1ivs', '1j1u', '1jbs', '1jbt', '1knz', '1kog', '1kq2', '1kuq', '1lng', '1m8v', '1m8w',
'1m8x', '1m8y', '1mms', '1n1h', '1n35', '1nb7', '1o0b', '1o0c', '1ooa', '1p6v', '1pgl', '1qa6', '1qrs', '1qrt', '1qru',
'1qtq', '1rc7', '1rpu', '1s03', '1sds', '1sj3', '1sj4', '1sjf', '1tfw', '1tfy', '1u0b', '1un6', '1uvi', '1uvj', '1uvl', '1uvm', '1uvn',
'1vbx', '1vby', '1vbz', '1vc0', '1vc6', '1vfg', '1wmq', '1wpu', '1wrq', '1xnq', '1xnr', '1yyk', '1yyo', '1yyw', '1yz9', '1zbh',
'1zdh', '1zdi', '1zdj', '1zdk', '1zjw', '1zse', '2ake', '2ann', '2asb', '2atw', '2b2e', '2b2g', '2bbv', '2bgg', '2bny', '2bq5',
'2bs0', '2bs1', '2bx2', '2c0b', '2c4r', '2c50', '2c51', '2d6f', '2db3', '2der', '2det', '2dr2', '2dr5', '2dr7', '2dr8', '2dr9',
'2dra', '2drb', '2du3', '2du4', '2du5', '2du6', '2dvi', '2e5l', '2e9r', '2e9z', '2ec0', '2ez6', '2f4v', '2f8k', '2f8s', '2f8t',
'2g4b', '2gic', '2hvy', '2hyi', '2i91', '2ix1', '2iy5', '2iz8', '2izm', '2izn', '2j0s', '2jlu', '2jlv', '2jlw', '2jlx', '2jly', '2jlz', '2nue',
'2nuf', '2nug', '2oih', '2oj3', '2ozb', '2pjp', '2ply', '2q66', '2qqp', '2qux', '2r7r', '2r7t', '2r7u', '2r7v', '2r7w', '2r7x', '2r8s',
'2r92', '2rd2', '2re8', '2rfk', '2uwm', '2uxb', '2v3c', '2vod', '2von', '2vop', '2w2h', '2xgj', '2xnr', '2xs2', '2xs5', '2xs7',
'2xxa', '2xzl', '2xzo', '2yhm', '2yjy', '2ykg', '2zh1', '2zh2', '2zh3', '2zh4', '2zh5', '2zh6', '2zh7', '2zh8', '2zh9', '2zha',
'2zhb', '2zko', '2zm6', '2zni', '2zxu', '2zzn', '3a2k', '3adi', '3adl', '3aev', '3ahu', '3akz', '3al0', '3am1', '3amu', '3avt',
'3avu', '3avv', '3avw', '3avx', '3avy', '3boy', '3bsb', '3bso', '3bsx', '3bt7', '3bx2', '3bx3', '3d2s', '3eph', '3epj', '3eqt',
'3ex7', '3fht', '3foz', '3ftf', '3g0h', '3g9y', '3gib', '3h5x', '3h5y', '3hhz', '3hjy', '3hsb', '3htx', '3i5x', '3i61', '3i62', '3iev',
'3ivk', '3iwn', '3jb7', '3k0j', '3k49', '3k4e', '3k5q', '3k5y', '3k5z', '3k61', '3k62', '3k64', '3klv', '3kmq', '3kms', '3kna',
'3koa', '3ks8', '3ktw', '3l25', '3l26', '3lwq', '3lwr', '3lwv', '3m7n', '3m85', '3mdg', '3mdi', '3mj0', '3mqk', '3ndb', '3nl0',
'3nmu', '3nnc', '3nnh', '3nvi', '3o3i', '3og8', '3oij', '3ouy', '3ov7', '3ova', '3ovb', '3ovs', '3pew', '3pey', '3pf4', '3pkm',
'3pla', '3ptx', '3pu0', '3pu1', '3pu4', '3q0l', '3q0m', '3q0n', '3q0o', '3q0p', '3q0q', '3q0r', '3q0s', '3q1q', '3q1r', '3q2t',
'3qg9', '3qgb', '3qgc', '3qjj', '3qjl', '3qjp', '3qsu', '3qsy', '3r2c', '3r2d', '3r9w', '3r9x', '3rc8', '3rer', '3siu', '3siv', '3sqw',
'3sqx', '3t3n', '3t3o', '3t5n', '3t5q', '3tup', '3u4m', '3u56', '3umy', '3v6y', '3v71', '3v74', '3v7e', '3vjr', '3vnu', '3vnv',
'3vyx', '3vyy', '3wbm', '3wc2', '3wqy', '3wqz', '3wzi', '3zd6', '3zd7', '3zla', '4a2x', '4al5', '4al6', '4al7', '4alp', '4b3g',
'4ba2', '4bhh', '4bpb', '4bw0', '4c4w', '4c7o', '4d25', '4d26', '4db2', '4ed5', '4f02', '4f1n', '4fsj', '4ftb', '4fte', '4fts',
'4fvu', '4fwt', '4g9z', '4gha', '4ghl', '4gl2', '4h5o', '4h5p', '4ht8', '4ig8', '4ilm', '4j1g', '4j39', '4j5v', '4j7l', '4jgn', '4jk0',
'4jng', '4jvy', '4jya', '4k4t', '4k50', '4khp', '4kji', '4knq', '4kq0', '4kr6', '4kr9', '4ktg', '4l8h', '4lg2', '4lmz', '4m2z', '4m30',
'4m4o', '4m6d', '4m7a', '4m7d', '4n2q', '4n2s', '4n48', '4ngb', '4ngc', '4ngd', '4ngf', '4ngg', '4nh3', '4nh5', '4nh6',
'4nl3', '4o8j', '4ohz', '4oo1', '4oog', '4ox9', '4p3e', '4pdb', '4pmi', '4pmw', '4prf', '4qg3', '4qi2', '4qil', '4qm6', '4qoz',
'4qpx', '4qqb', '4qvc', '4qvd', '4r3i', '4rcj', '4rwn', '4rwo', '4rwp', '4s3n', '4tuw', '4tux', '4tyw', '4tyy', '4tz0', '4tz6',
'4uft', '4w5q', '4wan', '4wc2', '4wc3', '4wc4', '4wc5', '4wc6', '4wc7', '4wj4', '4wkr', '4wrt', '4wsa', '4wsb', '4wta',
'4wtc', '4wtd', '4wte', '4wtf', '4wtg', '4wzm', '4wzq', '4x0a', '4x0b', '4x2b', '4x4o', '4xbf', '4xco', '4xww', '4y91',
'4yb1', '4yco', '4yhh', '4yvi', '4yvj', '4yvk', '4yye', '4z31', '4z7l', '4z92', '4zld', '4zlr', '4zt9', '5a0t', '5a0v', '5amq', '5amr',
'5aor', '5bud', '5bym', '5bz1', '5bz5', '5bzv', '5c9h', '5ccb', '5d6g', '5dar', '5ddo', '5ddp', '5de5', '5det', '5dno', '5dto',
'5dv7', '5e02', '5e08', '5e3h', '5eim', '5elh', '5elk', '5elr', '5els', '5elt', '5elx', '5emo', '5en1', '5ex7', '5f5f', '5f5h', '5f8g', '5f8h', '5f8i', '5f8j', '5f8l', '5f8m', '5f8n', '5f98', '5f9f', '5f9h', '5fj4', '5fmz', '5fn1', '5fvc', '5g4u', '5g4v', '5gin', '5gio',
'5gip', '5gjb', '5gxh', '5gxi', '5h0r', '5h1k', '5h1l', '5h3u', '5hab', '5hc9', '5ho4', '5hsw', '5i9d', '5i9f', '5i9g', '5i9h', '5iwa',
'5jaj', '5jbj', '5jc3', '5jc7', '5jcf', '5jch', '5jrc', '5jxs', '5k78', '5kl1', '5kla', '5l2l', '5lm7', '5lmn', '5lta', '5m0i', '5m0j',
'5m3h', '5m3j', '5mfx', '5mjv', '5msf', '5n94', '5npm', '5o1y', '5o1z', '5o20', '5o5j', '5oc6', '5sup', '5t8y', '5the', '5tsn',
'5ud5', '5udz', '5uj2', '5uz9', '5v6x', '5v7c', '5vw1', '5vzj', '5w0m', '5w0o', '5w1h', '5wea', '5wlh', '5ws2', '5wt3',
'5wty', '5wwe', '5wwf', '5wwg', '5www', '5wwx', '5wzg', '5wzh', '5wzi', '5wzj', '5wzk', '5xc6', '5xlo', '5xlp', '5y58',
'5y7m', '5yki', '5yts', '5ytt', '5ytv', '5ytx', '5z4d', '5z4j', '5z98', '5z9w', '5z9x', '5zc9', '5zw4', '6a6j', '6a6l', '6aax',
'6agb', '6ah3', '6ahr', '6ahu', '6b0b', '6b14', '6b3k', '6b45', '6b46', '6b47', '6bbo', '6bjg', '6bjh', '6bjv', '6bjy', '6cf2',
'6cmn', '6db8', '6db9', '6dcl', '6dnh', '6dtd', '6dti', '6du4', '6dzk', '6e9e', '6e9f', '6een', '6evk', '6f4a', '6f4g', '6f4h',
'6fhi', '6fpq', '6fpx', '6fq3', '6fql', '6fuw', '6g19', '6g1s', '6g1x', '6gd2', '6gd3', '6gjz', '6gkh', '6gkm', '6gpg', '6gv4',
'6gvy', '6gx6', '6h25', '6h5q', '6h5s', '6h61', '6h66', '6hct', '6htu', '6i0u', '6i0v', '6i3p', '6ifk', '6ifl', '6ifn', '6ifo', '6ifr', '6ifu', '6ify', '6ifz', '6ig0', '6iqw', '6iv6', '6iv9', '6jim', '6jvx', '6k0b', '6kl9', '6klb', '6kle', '6klh', '6ktc', '6kug', '6kuj', '6kuk', '6kup', '6kur', '6kut', '6kuu', '6kuv', '6kwq', '6kyv', '6l1w', '6llb', '6lnc', '6lse', '6lsf', '6lsg', '6lsh', '6m7d', '6mfn', '6mkn', '6msf', '6mur', '6mus', '6mut', '6muu', '6nm9', '6nma', '6nmc', '6nmd', '6noc', '6nod', '6nof', '6noh', '6nud', '6nue', '6nut', '6ny5', '6o16', '6o1k', '6o1l', '6o1m', '6o5f', '6o6v', '6o6x', '6o7b', '6o7k', '6p7m', '6pif', '6pig', '6ppn', '6ppp', '6ppq', '6ppv', '6pun', '6pzq', '6q8u', '6qic', '6qwl', '6qx3', '6qx8', '6qxe', '6r7g', '6r9j', '6r9o', '6r9p', '6r9q', '6ra4', '6scf', '6sjd', '6spc', '6spe', '6sx0', '6sx2', '6sy6', '6szu', '6szv', '6t0n', '6tqb', '6tz2', '6u8d', '6u8k', '6u9x', '6uro', '6uv1', '6uv2', '6uv3', '6uv4', '6v4x', '6v5b', '6v9q', '6vqv', '6vqw', '6vqx', '6vrb', '6vrc', '6w11', '6w1x', '6w62', '6w6k', '6wxx', '6wxy', '6x5m', '6x5n', '6xez', '6xh1', '6xh2', '6xh3', '6xki', '6xl1', '6xqb', '6yud', '6yym', '6yyt', '6z6b', '6z6g', '6z8k', '6zdq', '6zlc', '6zot', '6zww', '7aap', '7aav', '7abh', '7b3d', '7bg6', '7bv2', '7bzf', '7c45', '7cgf', '7cgg', '7cgh', '7cgi', '7cgj', '7cgk', '7cgl', '7cgm', '7ctt', '7cxm', '7cxn', '7cyq', '7d7v', '7dd3', '7dic', '7did', '7dmq', '7dol', '7el9', '7ela', '7elc', '7ewq', '7jhy', '7jl0', '7jl1', '7jl2', '7jl3', '7jzw', '7jzx', '7jzy', '7jzz', '7k98', '7k9m', '7ka0', '7kab', '7kha', '7krp', '7kwg', '7kx9', '7ld5', '7m5o', '7mlx', '7msf', '7n0b', '7n0c', '7n0d', '7nq4', '7nun', '7nuq', '7oe1', '7og6', '7ogk', '7oi3', '7om6', '7om7', '7oma', '7onb', '7onu', '7pmq']

a2021_r = ['7e9i', '7eem', '7efg', '7efh', '7eog', '7eoi', '7eok', '7eol', '7eom', '7eon', '7eoo', '7eop', '7mkt', '7ptq', '7qsh', '7sam',
'7sxp', '7tzt', '7uri', '7urm', '7vft', '7xd3', '7xd7', '7xsk', '7xsl', '7xsm', '7xsn', '7y2b', '8feo', '8fep', '8hb8', '8hd6',
'8hd7', '8i7n', '8keb', '8ked', '8q6n', '8rui', '8ruj', '8rul', '8rum', '8run', '8tkk', '8u5j', '8upt', '8upy', '8xtq', '9c6i', '9c6k',
'9cbu', '9cbw', '9cbx', '9cby', '9cxf', '9de8', '9ely', '9g7c']

b2021_r = ['1osu', '1q9a', '1z43', '255d', '283d', '2a43', '2eet', '2eev', '2g9c', '2gis', '2q1r', '2xnw', '2xo0', '2xo1', '2ygh', '2zy6',
'387d', '3a3a', '3d0m', '3dvz', '3fo4', '3gao', '3gca', '3gog', '3gx2', '3gx5', '3iqn', '3iqp', '3l0u', '3la5', '3nd3', '3nd4',
'3q50', '3r4f', '3sd3', '413d', '483d', '4aob', '4b5r', '4c40', '4cs1', '4ds6', '4e8m', '4e8n', '4e8p', '4e8q', '4e8v', '4en5',
'4faq', '4fe5', '4fej', '4fen', '4feo', '4fep', '4jiy', '4k27', '4l81', '4lvv', '4lvw', '4lvy', '4lvz', '4lw0', '4lx5', '4lx6', '4p8z',
'4rbq', '4tzx', '4tzy', '4xnr', '5c5w', '5da6', '5fk2', '5fk3', '5fk5', '5fke', '5fkf', '5g2y', '5m0h', '6d3p', '6e80', '6jxm',
'6pq7', '6t3s', '6tb7', '6tf1', '6tf3', '6tfe', '6ubu', '6uc7', '6uc8', '6ues', '6uet', '6vui', '6wjs', '6y0y', '7d7w', '7d81',
'7d82', '7eaf', '7kd1', '7lyj', '7mky']

a2021_rr = ['7q7x', '7q7y', '7q80', '7q82', '7qtn', '7xd5', '7xd6', '7y2p', '7ycg', '7ych', '7yci', '7yg8', '7yg9', '7yga', '7ygd', '8clm',
'8hb1', '8hb3', '8hzd', '8i3z', '8oly', '8olz', '8x9l', '8x9n', '8xtp', '8xtr', '8za0', '8za4', '9is7']

b2021_rr = ['1csl', '1d4r', '1dqf', '1dqh', '1duq', '1ik5', '1jzv', '1lnt', '1mhk', '1nta', '1ntb', '1nuj', '1nuv', '1nyi', '1q29', '1qc0', '1qcu', '1r3o', '1sdr', '1t0d', '1t0e', '1x9c', '1x9k', '1ykq', '1ykv', '1yls', '1zft', '1zfv', '1zfx', '1zx7', '1zz5', '2a04', '2bcy', '2bcz', '2d2k', '2f4t', '2fgp', '2g32', '2gcs', '2gcv', '2gpm', '2gq4', '2gq5', '2gq6', '2gq7', '2h0w', '2h0x', '2ho6', '2ho7', '2jlt', '2nok', '2npy', '2npz', '2oeu', '2oue', '2pn3', '2pn4', '2quw', '2v6w', '2v7r', '2val', '2vuq', '2w89', '2x2q', '2xsl', '2yie',
'2yif', '353d', '354d', '357d', '364d', '397d', '3b31', '3b4a', '3b4b', '3b4c', '3bnp', '3bnq', '3bnr', '3bns', '3cgp', '3cgq',
'3cgr', '3cgs', '3eog', '3f4e', '3f4g', '3f4h', '3fs0', '3ftm', '3g78', '3gs1', '3gs8', '3gvn', '3igi', '3loa', '3mj3', '3mja',
'3mjb', '3owi', '3ox0', '3oxe', '3oxm', '3p22', '3sj2', '3ski', '3tzr', '3zd4', '3zd5', '3zp8', '429d', '434d', '435d', '439d',
'464d', '466d', '472d', '4e6b', '4e8k', '4e8t', '4f8u', '4far', '4faw', '4jrt', '4kz2', '4mce', '4mcf', '4meg', '4meh', '4mgn',
'4p97', '4phy', '4qjd', '4qjh', '4r0d', '4u34', '4u35', '4u37', '4u38', '5c45', '5cnr', '5dh7', '5g4t', '5j02', '5ktj', '5kvj',
'5kx9', '5t3k', '5tko', '5uef', '5ux3', '5v0o', '5v9z', '5vcf', '5vci', '5vgw', '5zei', '6az4', '6bfb', '6c63', '6c64', '6c65',
'6c8m', '6chr', '6cih', '6dn1', '6dn2', '6dn3', '6hu6', '6jbg', '6owl', '6r47', '6t3k', '6t3r', '6u7z', '6u8f', '6ufh', '6ufj', '6ufk',
'6ufm', '6yml', '6ymm', '6zpf', '6zq9', '6zr1', '6zrl', '6zrs', '6zw3', '6zwu', '6zx5', '6zx8', '7a9l', '7a9n', '7a9o', '7a9p',
'7a9q', '7a9r', '7a9s', '7a9t', '7ez2']

valid_ids = a2021_rr + b2021_rr + a2021_r + b2021_r + a2021_rp + b2021_rp
rna_rna = a2021_rr + b2021_rr + a2021_r + b2021_r
protein_rna = a2021_rp + b2021_rp

def download_pdb(pdb_id, save_dir="pdb_structures", file_format="cif"):
    os.makedirs(save_dir, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
    dst = Path(save_dir) / f"{pdb_id}_gt.{file_format}"
    print(f"⬇️ Downloading GT structure {pdb_id} to {dst}")
    urllib.request.urlretrieve(url, dst)
    print(f"✅ Download complete: {dst}")
    return dst

parser = argparse.ArgumentParser(description="Prepare evaluation data for RNA-RNA and Protein-RNA predictions.")
parser.add_argument('--pred_data_dir', type=str, default='evaluation/predictions/rnaformer',
                    help='Directory containing prediction data in CIF format.')
parser.add_argument('--gt_data_dir', type=str, default='evaluation/predictions/gt',
                    help='Directory containing ground truth data in CIF format.')
args = parser.parse_args()

# pred_data_dir = Path('/home/fred/current_projects/github/RnaBench/synthetic_msa_nature_methods/af3_complex_cifs')

pred_data_dir = Path(args.pred_data_dir)

gt_data_dir = Path(args.gt_data_dir)

# rna_rna_pred_outdir = Path('rna_rna', 'af')
# protein_rna_pred_outdir = Path('protein_rna', 'af')

rna_rna_pred_outdir = Path('data', 'rna_rna', 'af')
protein_rna_pred_outdir = Path('data', 'protein_rna', 'af')

rna_rna_pred_outdir.mkdir(exist_ok=True, parents=True)
protein_rna_pred_outdir.mkdir(exist_ok=True, parents=True)

rna_rna_gt_outdir = Path('data', 'rna_rna', 'pdb')
protein_rna_gt_outdir = Path('data', 'protein_rna', 'pdb')

rna_rna_gt_outdir.mkdir(exist_ok=True, parents=True)
protein_rna_gt_outdir.mkdir(exist_ok=True, parents=True)

raw_pred_data = list(pred_data_dir.glob('*.cif'))
raw_gt_data = list(gt_data_dir.glob('*.cif'))

print(len(raw_pred_data))
print(len(raw_gt_data))

raw_pred_ids = [str(p.stem)[:4].lower() for p in raw_pred_data]
raw_gt_ids = [str(p.stem)[:4].lower() for p in raw_gt_data]

for i in raw_pred_ids:
    if not i in raw_gt_ids:
        download_pdb(i, save_dir=gt_data_dir, file_format="cif")

raw_gt_data = list(gt_data_dir.glob('*.cif'))

pred_data = []
gt_data = []

rna_rna_pred_data = []
protein_rna_pred_data = []

rna_rna_gt_data = []
protein_rna_gt_data = []

for p in raw_pred_data:
    # pdb_id = str(p.stem).split('_')[1].lower()
    pdb_id = str(p.stem)[:4].lower()
    if not pdb_id.lower() in valid_ids:
        print('Invalid data point:', pdb_id)
    else:
        if pdb_id in rna_rna:
            rna_rna_pred_data.append(p)
            dest = Path(rna_rna_pred_outdir, f"fold_{pdb_id}_s1_model_0.cif")
            shutil.copy(p, dest)
        else:
            protein_rna_pred_data.append(p)
            dest = Path(protein_rna_pred_outdir, f"fold_{pdb_id}_s1_model_0.cif")
            shutil.copy(p, dest)

        pred_data.append(p)

for p in raw_gt_data:
    pdb_id = str(p.stem)[:4].lower()
    if not pdb_id in valid_ids:
        print('Invalid data point:', pdb_id)
    else:
        if pdb_id in rna_rna:
            rna_rna_gt_data.append(p)
            dest = Path(rna_rna_gt_outdir, f"{pdb_id}.cif")
            shutil.copy(p, dest)
        else:
            protein_rna_gt_data.append(p)
            dest = Path(protein_rna_gt_outdir, f"{pdb_id}.cif")
            shutil.copy(p, dest)
        gt_data.append(p)

print(len(pred_data))
print(len(gt_data))

print(len(rna_rna_pred_data))
print(len(rna_rna_gt_data))
print(len(protein_rna_pred_data))
print(len(protein_rna_gt_data))