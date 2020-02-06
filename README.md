# FormulaeforElasticityinNestedLogitModel
{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Formulae Auto-Drivation for Nested Logit Models\n",
    "\n",
    "This file uses the Symbolic Python (SymPy) package to derive the elasticity formulae for Nested Logit Models. SymPy is a Python library for symbolic mathematics. It aims to become a full-featured computer algebra system (CAS) while keeping the code as simple as possible in order to be comprehensible and easily extensible. More details about SymPy can be found on the SymPy package main page: [SymPy](https://www.sympy.org/en/index.html).\n",
    "\n",
    "To use the SymPy package, we first have to import the SymPy package as follows. Note that we also include the `init_printing()` function to enable the pretty print which is particularly useful when we are printing a large amount of mathematical formulae."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "#### Import libraries\n",
    "from sympy import *\n",
    "import sys\n",
    "init_printing()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Next, let's set up our model by specifying a few basic parameters in the model."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "#### Settings\n",
    "number_end_group = 3  # This is the number of vehicles in each end group of our tree\n",
    "number_splits_non_end_node = 2  # This is the number of splits for each node if the node is not on the last two levels\n",
    "depth_of_tree = 3  # This is the depth of our tree"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "From the above settings, we build up a tree with the following structre."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 49,
   "metadata": {
    "tags": [
     "hide_input"
    ]
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "%matplotlib inline\n",
    "\n",
    "decision_node = dict(boxstyle=\"sawtooth\", fc=\"0.8\")\n",
    "leaf_node = dict(boxstyle=\"round4\", fc=\"0.8\")\n",
    "arrow_args = dict(arrowstyle=\"<-\")\n",
    "\n",
    "\n",
    "def get_leaf_num(tree):\n",
    "    leaf_num = 0\n",
    "    for key in tree.keys():\n",
    "        if type(tree[key]).__name__ == \"dict\":\n",
    "            leaf_num +=get_leaf_num(tree[key])\n",
    "        else:\n",
    "            leaf_num +=1\n",
    "    return leaf_num\n",
    "\n",
    "\n",
    "def get_tree_depth(tree):\n",
    "    depth = 0\n",
    "    for key in tree.keys():\n",
    "        if type(tree[key]).__name__ == \"dict\":\n",
    "            thisdepth = 1+ get_tree_depth(tree[key])\n",
    "        else:\n",
    "            thisdepth = 1\n",
    "        if thisdepth>depth: depth = thisdepth\n",
    "    return depth\n",
    "\n",
    "\n",
    "def plotNode(nodeTxt, centerPt, parentPt, nodeType):\n",
    "    createPlot.ax1.annotate(nodeTxt, xy=parentPt, xycoords='axes fraction',\n",
    "                            xytext=centerPt, textcoords='axes fraction',\n",
    "                            va=\"center\", ha=\"center\", bbox=nodeType, arrowprops=arrow_args)\n",
    "\n",
    "    \n",
    "def plotTree(myTree, parentPt, nodeTxt):\n",
    "    numLeafs = get_leaf_num(myTree)\n",
    "    depth = get_tree_depth(myTree)\n",
    "    firstStr = list(myTree.keys())[0]\n",
    "    cntrPt = (plotTree.xOff + (1.0 + float(numLeafs)) / 2.0 / plotTree.totalW, plotTree.yOff)\n",
    "    plotNode(firstStr, cntrPt, parentPt, decision_node)\n",
    "    secondDict = myTree[firstStr]\n",
    "    plotTree.yOff = plotTree.yOff - 1.0 / plotTree.totalD\n",
    "    for key in myTree.keys():\n",
    "        if type(myTree[key]).__name__ == 'dict':\n",
    "            plotTree(myTree[key], cntrPt, str(key))\n",
    "        else:\n",
    "            plotTree.xOff = plotTree.xOff + 1.0 / plotTree.totalW\n",
    "            plotNode(myTree[key], (plotTree.xOff, plotTree.yOff), cntrPt, leaf_node)\n",
    "    plotTree.yOff = plotTree.yOff + 1.0 / plotTree.totalD\n",
    "\n",
    "\n",
    "def createPlot(inTree):\n",
    "    fig = plt.figure(1, facecolor='white')\n",
    "    fig.clf()\n",
    "    axprops = dict(xticks=[], yticks=[])\n",
    "    createPlot.ax1 = plt.subplot(111, frameon=False, **axprops)\n",
    "    plotTree.totalW = float(get_leaf_num(inTree))\n",
    "    plotTree.totalD = float(get_tree_depth(inTree))\n",
    "    plotTree.xOff = -0.5 / plotTree.totalW\n",
    "    plotTree.yOff = 1.0\n",
    "    plotTree(inTree, (0.5, 1.0), '')\n",
    "    plt.show()\n",
    "\n",
    "\n",
    "def create_dict(level):\n",
    "    global count, number_end_group, number_splits_non_end_node, depth_of_tree\n",
    "    if level == depth_of_tree - 1:\n",
    "        dic = dict()\n",
    "        for i in range(number_end_group):\n",
    "            dic[str(i)] = count\n",
    "            count += 1\n",
    "        return dic\n",
    "    else:\n",
    "        dic = dict()\n",
    "        for i in range(number_splits_non_end_node):\n",
    "            name = str(i)\n",
    "            dic[str(i)] = create_dict(level+1)\n",
    "        return dic\n",
    "\n",
    "    \n",
    "count = 0\n",
    "createPlot(create_dict(0))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now that we finished setting up our model we can start calculating the formulae of interest. Let's first review the formulae for the share and elasticities in the Nested Logit Model. The share of one end node in the Nested Logit Model is\n",
    "$$ s_j = \\frac{\\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right)^{\\theta_{0 0}} \\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} \\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} e^{\\frac{v_{0}}{\\theta_{2 0}}}}{\\left(\\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right)^{\\theta_{0 0}} + 1\\right) \\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right) \\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right) \\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)}. $$ "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Self-elasticity\n",
    "The self-elasticity of the node is\n",
    "$$ self\\_elas = \\beta p_{0} \\left(- s_{3} + \\frac{1}{\\theta_{2 0}} - \\frac{s_{3}}{s_{2} \\theta_{2 0}} + \\frac{s_{3}}{s_{2} \\theta_{1 0}} - \\frac{s_{3}}{s_{1} \\theta_{1 0}} + \\frac{s_{3}}{s_{1} \\theta_{0 0}} + \\frac{s_{3}}{s_{0}} - \\frac{s_{3}}{s_{0} \\theta_{0 0}}\\right), $$ and the cross-elasticity of the node with another node within the same group in the second to last layer is\n",
    "$$ cross\\_elas = \\beta p_{1} \\left(- s_{3} - \\frac{s_{3}}{s_{2} \\theta_{2 0}} + \\frac{s_{3}}{s_{2} \\theta_{1 0}} - \\frac{s_{3}}{s_{1} \\theta_{1 0}} + \\frac{s_{3}}{s_{1} \\theta_{0 0}} + \\frac{s_{3}}{s_{0}} - \\frac{s_{3}}{s_{0} \\theta_{0 0}}\\right). $$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Basic elements for calculation\n",
    "total_end_nodes = number_splits_non_end_node ** (depth_of_tree - 1) * number_end_group\n",
    "for i in range(total_end_nodes):\n",
    "    exec('p_%d = symbols(\"p_%d\")' % (i, i))\n",
    "total_theta = 0\n",
    "for i in range(depth_of_tree):\n",
    "    total_theta += number_splits_non_end_node ** i\n",
    "for i in range(total_end_nodes):\n",
    "    exec('v_%d = symbols(\"v_%d\")' % (i, i))\n",
    "for i in range(depth_of_tree):\n",
    "    for j in range(number_splits_non_end_node ** i):\n",
    "        exec('theta_%d_%d = symbols(\"theta_%d_%d\")' % (i, j, i, j))\n",
    "for i in range(depth_of_tree+1):\n",
    "    exec('s_%d = symbols(\"s_%d\")' % (i, i))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Intermediate sums in s\n",
    "for i in range(total_end_nodes):\n",
    "    exec('sum_exp_%d_%d = exp(v_%d / theta_%d_%d)' % (depth_of_tree, i, i, depth_of_tree-1, i // number_end_group))\n",
    "    \n",
    "\n",
    "for i in range(total_end_nodes // number_end_group):\n",
    "    st = 'sum_%d_%d =' % (depth_of_tree-1, i)\n",
    "    for j in range(i*number_end_group, i*number_end_group + number_end_group):\n",
    "        st = st + ' exp(v_%d / theta_%d_%d) +' % (j, depth_of_tree-1, i)\n",
    "    exec(st[:-1])\n",
    "    exec('sum_exp_%d_%d = sum_%d_%d ** (theta_%d_%d / theta_%d_%d)' % (depth_of_tree-1, i, depth_of_tree-1, i, depth_of_tree-1, i,\n",
    "                                                                       depth_of_tree-2, i//number_splits_non_end_node))\n",
    "\n",
    "for k in range(depth_of_tree-2, 0, -1):\n",
    "    for i in range(number_splits_non_end_node ** k):\n",
    "        st = 'sum_%d_%d =' % (k, i)\n",
    "        for j in range(i*number_splits_non_end_node, (i+1)*number_splits_non_end_node):\n",
    "            st = st + ' sum_exp_%d_%d +' % (k+1, j)\n",
    "        exec(st[:-1])\n",
    "        exec('sum_exp_%d_%d = sum_%d_%d ** (theta_%d_%d / theta_%d_%d)' % (k, i, k, i, k, i,\n",
    "                                                                       k-1, i//number_splits_non_end_node))\n",
    "\n",
    "st = 'sum_0_0 ='\n",
    "for j in range(number_splits_non_end_node):\n",
    "    st = st + ' sum_exp_%d_%d +' % (1, j)\n",
    "exec(st[:-1])\n",
    "exec('sum_exp_0_0 = sum_0_0 ** (theta_0_0)')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Calculate final expression for s\n",
    "exec('s = sum_exp_%d_0 / (1 + sum_0_0 ** theta_0_0)' % depth_of_tree)\n",
    "for i in range(depth_of_tree):\n",
    "    exec('s = s / sum_%d_0 * sum_exp_%d_0' % (i, i))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/latex": [
       "$$\\frac{\\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right)^{\\theta_{0 0}} \\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} \\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} e^{\\frac{v_{0}}{\\theta_{2 0}}}}{\\left(\\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right)^{\\theta_{0 0}} + 1\\right) \\left(\\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right)^{\\frac{\\theta_{1 0}}{\\theta_{0 0}}} + \\left(\\left(e^{\\frac{v_{6}}{\\theta_{2 2}}} + e^{\\frac{v_{7}}{\\theta_{2 2}}} + e^{\\frac{v_{8}}{\\theta_{2 2}}}\\right)^{\\frac{\\theta_{2 2}}{\\theta_{1 1}}} + \\left(e^{\\frac{v_{10}}{\\theta_{2 3}}} + e^{\\frac{v_{11}}{\\theta_{2 3}}} + e^{\\frac{v_{9}}{\\theta_{2 3}}}\\right)^{\\frac{\\theta_{2 3}}{\\theta_{1 1}}}\\right)^{\\frac{\\theta_{1 1}}{\\theta_{0 0}}}\\right) \\left(\\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)^{\\frac{\\theta_{2 0}}{\\theta_{1 0}}} + \\left(e^{\\frac{v_{3}}{\\theta_{2 1}}} + e^{\\frac{v_{4}}{\\theta_{2 1}}} + e^{\\frac{v_{5}}{\\theta_{2 1}}}\\right)^{\\frac{\\theta_{2 1}}{\\theta_{1 0}}}\\right) \\left(e^{\\frac{v_{0}}{\\theta_{2 0}}} + e^{\\frac{v_{1}}{\\theta_{2 0}}} + e^{\\frac{v_{2}}{\\theta_{2 0}}}\\right)}$$"
      ],
      "text/plain": [
       "                                                                              \n",
       "                                                              ⎛               \n",
       "                                                              ⎜               \n",
       "                                                              ⎜               \n",
       "                                                              ⎜⎛              \n",
       "                                                              ⎜⎜              \n",
       "                                                              ⎜⎜              \n",
       "                                                              ⎜⎜⎛  v₀      v₁ \n",
       "                                                              ⎜⎜⎜ ────    ────\n",
       "                                                              ⎜⎜⎜ θ₂ ₀    θ₂ ₀\n",
       "                                                              ⎝⎝⎝ℯ     + ℯ    \n",
       "──────────────────────────────────────────────────────────────────────────────\n",
       "⎛                                                                             \n",
       "⎜⎛                                                           θ₁ ₀             \n",
       "⎜⎜                                                           ────             \n",
       "⎜⎜                                                           θ₀ ₀             \n",
       "⎜⎜⎛                       θ₂ ₀                          θ₂ ₁⎞       ⎛         \n",
       "⎜⎜⎜                       ────                          ────⎟       ⎜         \n",
       "⎜⎜⎜                       θ₁ ₀                          θ₁ ₀⎟       ⎜         \n",
       "⎜⎜⎜⎛  v₀      v₁      v₂ ⎞       ⎛  v₃      v₄      v₅ ⎞    ⎟       ⎜⎛  v₆    \n",
       "⎜⎜⎜⎜ ────    ────    ────⎟       ⎜ ────    ────    ────⎟    ⎟       ⎜⎜ ────   \n",
       "⎜⎜⎜⎜ θ₂ ₀    θ₂ ₀    θ₂ ₀⎟       ⎜ θ₂ ₁    θ₂ ₁    θ₂ ₁⎟    ⎟       ⎜⎜ θ₂ ₂   \n",
       "⎝⎝⎝⎝ℯ     + ℯ     + ℯ    ⎠     + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠     + ⎝⎝ℯ     + \n",
       "\n",
       "                                                                              \n",
       "                                            θ₁ ₀                              \n",
       "                                            ────                              \n",
       "                                            θ₀ ₀                              \n",
       "         θ₂ ₀                          θ₂ ₁⎞       ⎛                       θ₂ \n",
       "         ────                          ────⎟       ⎜                       ───\n",
       "         θ₁ ₀                          θ₁ ₀⎟       ⎜                       θ₁ \n",
       "     v₂ ⎞       ⎛  v₃      v₄      v₅ ⎞    ⎟       ⎜⎛  v₆      v₇      v₈ ⎞   \n",
       "    ────⎟       ⎜ ────    ────    ────⎟    ⎟       ⎜⎜ ────    ────    ────⎟   \n",
       "    θ₂ ₀⎟       ⎜ θ₂ ₁    θ₂ ₁    θ₂ ₁⎟    ⎟       ⎜⎜ θ₂ ₂    θ₂ ₂    θ₂ ₂⎟   \n",
       " + ℯ    ⎠     + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠     + ⎝⎝ℯ     + ℯ     + ℯ    ⎠   \n",
       "──────────────────────────────────────────────────────────────────────────────\n",
       "                                                      θ₀ ₀    ⎞               \n",
       "                                                 θ₁ ₁⎞        ⎟ ⎛             \n",
       "                                                 ────⎟        ⎟ ⎜             \n",
       "                                                 θ₀ ₀⎟        ⎟ ⎜             \n",
       "              θ₂ ₂                          θ₂ ₃⎞    ⎟        ⎟ ⎜⎛            \n",
       "              ────                          ────⎟    ⎟        ⎟ ⎜⎜            \n",
       "              θ₁ ₁                          θ₁ ₁⎟    ⎟        ⎟ ⎜⎜            \n",
       "  v₇      v₈ ⎞       ⎛ v₁₀     v₁₁      v₉ ⎞    ⎟    ⎟        ⎟ ⎜⎜⎛  v₀      v\n",
       " ────    ────⎟       ⎜ ────    ────    ────⎟    ⎟    ⎟        ⎟ ⎜⎜⎜ ────    ──\n",
       " θ₂ ₂    θ₂ ₂⎟       ⎜ θ₂ ₃    θ₂ ₃    θ₂ ₃⎟    ⎟    ⎟        ⎟ ⎜⎜⎜ θ₂ ₀    θ₂\n",
       "ℯ     + ℯ    ⎠     + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠    ⎠     + 1⎠⋅⎝⎝⎝ℯ     + ℯ  \n",
       "\n",
       "                                     θ₀ ₀                                     \n",
       "                                θ₁ ₁⎞                                         \n",
       "                                ────⎟                                         \n",
       "                                θ₀ ₀⎟                                         \n",
       "₂                          θ₂ ₃⎞    ⎟     ⎛                       θ₂ ₀        \n",
       "─                          ────⎟    ⎟     ⎜                       ────        \n",
       "₁                          θ₁ ₁⎟    ⎟     ⎜                       θ₁ ₀        \n",
       "    ⎛ v₁₀     v₁₁      v₉ ⎞    ⎟    ⎟     ⎜⎛  v₀      v₁      v₂ ⎞       ⎛  v₃\n",
       "    ⎜ ────    ────    ────⎟    ⎟    ⎟     ⎜⎜ ────    ────    ────⎟       ⎜ ───\n",
       "    ⎜ θ₂ ₃    θ₂ ₃    θ₂ ₃⎟    ⎟    ⎟     ⎜⎜ θ₂ ₀    θ₂ ₀    θ₂ ₀⎟       ⎜ θ₂ \n",
       "  + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠    ⎠    ⋅⎝⎝ℯ     + ℯ     + ℯ    ⎠     + ⎝ℯ   \n",
       "──────────────────────────────────────────────────────────────────────────────\n",
       "                                                                              \n",
       "                                              θ₁ ₀                            \n",
       "                                              ────                            \n",
       "                                              θ₀ ₀                            \n",
       "           θ₂ ₀                          θ₂ ₁⎞       ⎛                       θ\n",
       "           ────                          ────⎟       ⎜                       ─\n",
       "           θ₁ ₀                          θ₁ ₀⎟       ⎜                       θ\n",
       "₁      v₂ ⎞       ⎛  v₃      v₄      v₅ ⎞    ⎟       ⎜⎛  v₆      v₇      v₈ ⎞ \n",
       "──    ────⎟       ⎜ ────    ────    ────⎟    ⎟       ⎜⎜ ────    ────    ────⎟ \n",
       " ₀    θ₂ ₀⎟       ⎜ θ₂ ₁    θ₂ ₁    θ₂ ₁⎟    ⎟       ⎜⎜ θ₂ ₂    θ₂ ₂    θ₂ ₂⎟ \n",
       "   + ℯ    ⎠     + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠     + ⎝⎝ℯ     + ℯ     + ℯ    ⎠ \n",
       "\n",
       "                                                                              \n",
       "                       θ₁ ₀                                                   \n",
       "                       ────                                                   \n",
       "                       θ₀ ₀                                                   \n",
       "                  θ₂ ₁⎞                            θ₂ ₀                       \n",
       "                  ────⎟                            ────                       \n",
       "                  θ₁ ₀⎟                            θ₁ ₀                       \n",
       "      v₄      v₅ ⎞    ⎟     ⎛  v₀      v₁      v₂ ⎞       v₀                  \n",
       "─    ────    ────⎟    ⎟     ⎜ ────    ────    ────⎟      ────                 \n",
       "₁    θ₂ ₁    θ₂ ₁⎟    ⎟     ⎜ θ₂ ₀    θ₂ ₀    θ₂ ₀⎟      θ₂ ₀                 \n",
       "  + ℯ     + ℯ    ⎠    ⎠    ⋅⎝ℯ     + ℯ     + ℯ    ⎠    ⋅ℯ                     \n",
       "──────────────────────────────────────────────────────────────────────────────\n",
       "                                                                              \n",
       "                                  θ₁ ₁⎞                                       \n",
       "                                  ────⎟                                       \n",
       "                                  θ₀ ₀⎟                                       \n",
       "₂ ₂                          θ₂ ₃⎞    ⎟ ⎛                       θ₂ ₀          \n",
       "───                          ────⎟    ⎟ ⎜                       ────          \n",
       "₁ ₁                          θ₁ ₁⎟    ⎟ ⎜                       θ₁ ₀          \n",
       "      ⎛ v₁₀     v₁₁      v₉ ⎞    ⎟    ⎟ ⎜⎛  v₀      v₁      v₂ ⎞       ⎛  v₃  \n",
       "      ⎜ ────    ────    ────⎟    ⎟    ⎟ ⎜⎜ ────    ────    ────⎟       ⎜ ──── \n",
       "      ⎜ θ₂ ₃    θ₂ ₃    θ₂ ₃⎟    ⎟    ⎟ ⎜⎜ θ₂ ₀    θ₂ ₀    θ₂ ₀⎟       ⎜ θ₂ ₁ \n",
       "    + ⎝ℯ     + ℯ     + ℯ    ⎠    ⎠    ⎠⋅⎝⎝ℯ     + ℯ     + ℯ    ⎠     + ⎝ℯ     \n",
       "\n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "─────────────────────────────────────────────\n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                                             \n",
       "                θ₂ ₁⎞                        \n",
       "                ────⎟                        \n",
       "                θ₁ ₀⎟                        \n",
       "    v₄      v₅ ⎞    ⎟ ⎛  v₀      v₁      v₂ ⎞\n",
       "   ────    ────⎟    ⎟ ⎜ ────    ────    ────⎟\n",
       "   θ₂ ₁    θ₂ ₁⎟    ⎟ ⎜ θ₂ ₀    θ₂ ₀    θ₂ ₀⎟\n",
       "+ ℯ     + ℯ    ⎠    ⎠⋅⎝ℯ     + ℯ     + ℯ    ⎠"
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "s"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [],
   "source": [
    "#### Deriving the own-elasticity formular\n",
    "s_diff = diff(s, v_0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAgsAAAAxBAMAAACmD/cfAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAEM3dMol2u0SZImZUq++CbGFdAAAACXBIWXMAAA7EAAAOxAGVKw4bAAAGuElEQVRoBb1aTWhcVRQ+Y2byM0kzg0VEBBMQunXCdCOiDDQW7CaxJeOqdFAzUEprhOImiwY3RdykFLpwNRvNShJFXLjQQkHQhVasGyVYsehGNIFaY/0Zz7lv5t337j3nnos85i0m953zne875+a+mTlzL0BB10Onv4pgKq3s1iNgcWQKStNSwiPSZCCbcJWxuqYpmGq4NuY+jkxBaVpKOJNWhOkAtls6rLZRvaujII5MQWlaSnhEmgzkvahpmK5HTUMcmYLStJRwpsYo06koVCXmoQCII1NQmpYSHlWPCyovuBb2/ljMWyTEkWkoRUsLZ/PXjFtR9VVf1XiMP45MQWlaSnhUoi6oerPcc23M/btwgrG6pjgyDaVoaeFuUlH3W4vPRSyH6sXFKxF0cWQKStNSwr0026tKgZXu8/Bbv+8FOgbkme339xyrextHpqEULS3cTQrvJ1qTm4w5Y7oIb2TupKHOYyLjyBSUpqWEcyVsw8RO9dz5FucztuoerItO6yAeeyeN4sg0lKKlhXPJba/34NDy+G3OZ2zVPz4SfRkH8ehXHJmGUrSUcLbbqPR/h3KrMi+XcKa/DO98eEsGGA/xBFdVEh9HpqA0rXA4321MXsK3yLHQgr68D9/ABWUaAHmCq2oQH0emoDStYDjXbeA30lodxo/IReJ30n1YhJMygjyGJ7yqCBZHpqA0LSWc6zYe6MBnmN7DHUqSvY7ADH2UPMU6U+OAJ7iqEBxHpqA0LSWca2nK3S49D2W5QT7ePYeAFzv4ErgSntCqMsFxZApK01LCxZbmsVb1r0CJ5Cq9pgASd2BVZeLjyDSUoiWHS93GWH36h0yW7HCpxZodY2BVZZFxZApK05LCbbfxZjYnfCTa53t5i3NXugPbzmeJQ2ECIlYV4nyySYebyHwU5CQVLSbcpIgvttv4dGiK/bsAH9fzWI4iZlUhi0c208lzmzsPBTlJTcsPH2rY/ihHOHSH/j7b/dZxcxTqqko4PDJ2GjxUfho0LT+cxKnpshdXg/VGjQqgSHXYaUi96aAAyXzTVQBhARRpfSObBqfpKqCGAihGNw0PNuk6Ctmm62qz+Wuz+TQl0c9daCg9bgJavjNF7qFPpMihMmS11JEfINdUs/nMj83mzf8jSTHDEnGY507vKOPBlW+6CvhXFkAxzA1G9lAA5JquAmoogGL00+A0XQXUUADF6KfBaboKqKEAitFPg9N0FVBDARSjn4ZEMW26bA3jp1tpNoHB6i3XaSmsx0dZnx35kuxbpE/GSVpaZ+SHW8Cw6bJNynqpYd3iaGLnsOuzFKmHQaW+zMCX5ForhoyRzNDmh0z4AMA1XZU1kH9wscQvQK1l76RRHKpQSSkVCOTCNF1LG/C3SJU6ygdQu57eSYM4FBQpKaUCoVyYputrKGu/O6HUoQYs7YiSQ0ccCoqUHEp7fyNzGcSV7sHMgceRGobbm3Md+HIjtbqDONSwxS1U0k0F4nJxw8bvQuW+a0zv063CV1pwppeanUEcCoYtbqGSTip2QzaYsRsFE//sfn/bHK1bWfWcQFuFxnlpd/dfoD2vSnvZg+VRpffxtw0GZVpcchpJeHQeQJI0TiNJRByZ2cMUw8kRDPdKmF2DuQ4drat0Jtc8L20VmnN3PwMuZdrzasNlBXX8JWBRpsUlp5EsdxsgSpITjCTJSZJiuHEEw70Sah34pE5H62Zb4D8ctFVozt3dg+l9c4zvDtR6LkkeBTcAOBRQi0tOIwkzDRAlyQlGkog4MpIUw40jGO5WAHMb8CTQ0brtHvzieWlbkpylAxjbBNzzOnEfZq97sBwKK62yKGpxaRqMJFUqSpLTSBLR5ywZSorh5FDC3RJqGzO0QVFpLNX93blke5N2Jw/gizrCTpX3/T1fF3UDOFSyeWlWg5HESmVJWg0kSUQfSJJiuHGEwt1JwC8EO4dbaD1Wx9ifXPdgqxDP3X1XOorO8gLmNdtxYC7KTIOHGmxe4mpIJJNpECRpGkiS5F6WJMWMjSMU7hSAt6Vrb+ErHq3DleR9m0y2Cunc3SNnewjbquMq9b5NuijzUHgoSFpcnIZEEiuVJWkaSJLk8KHwyIykGG4coXAshb3waJ15X5GciZ32vLj3K+O1p/OwUhFFD0VyYaWyJE2DuaS3SHKK4dYRCk8Esq90tA4/ZeaztnRsz93RnlcbXk892YFFUaUSyrxFJmFYqShpPikMjIgkMjHcOkLh2eyTsTnGt3K27nvQYs/40Z5XpX1FQb39xDURBeQ0V/Xknx2QJI3TwEhOkhTDrSMN/w9knkETJoUGEwAAAABJRU5ErkJggg==\n",
      "text/latex": [
       "$$- s_{3}^{2} + \\frac{s_{3}}{\\theta_{2 0}} - \\frac{s_{3}^{2}}{s_{2} \\theta_{2 0}} + \\frac{s_{3}^{2}}{s_{2} \\theta_{1 0}} - \\frac{s_{3}^{2}}{s_{1} \\theta_{1 0}} + \\frac{s_{3}^{2}}{s_{1} \\theta_{0 0}} + \\frac{s_{3}^{2}}{s_{0}} - \\frac{s_{3}^{2}}{s_{0} \\theta_{0 0}}$$"
      ],
      "text/plain": [
       "                   2         2         2         2       2       2  \n",
       "    2    s₃      s₃        s₃        s₃        s₃      s₃      s₃   \n",
       "- s₃  + ──── - ─────── + ─────── - ─────── + ─────── + ─── - ───────\n",
       "        θ₂ ₀   s₂⋅θ₂ ₀   s₂⋅θ₁ ₀   s₁⋅θ₁ ₀   s₁⋅θ₀ ₀    s₀   s₀⋅θ₀ ₀"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "s_simp = s_diff.subs(sum_exp_0_0 / (1+sum_exp_0_0), s_0)\n",
    "for i in range(depth_of_tree):\n",
    "    exec('s_simp = s_simp.subs(sum_exp_%d_0 / sum_%d_0, s_%d / s_%d)' % (i+1, i, i+1, i))\n",
    "s_simp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/latex": [
       "$$\\beta p_{0} \\left(- s_{3} + \\frac{1}{\\theta_{2 0}} - \\frac{s_{3}}{s_{2} \\theta_{2 0}} + \\frac{s_{3}}{s_{2} \\theta_{1 0}} - \\frac{s_{3}}{s_{1} \\theta_{1 0}} + \\frac{s_{3}}{s_{1} \\theta_{0 0}} + \\frac{s_{3}}{s_{0}} - \\frac{s_{3}}{s_{0} \\theta_{0 0}}\\right)$$"
      ],
      "text/plain": [
       "     ⎛       1        s₃        s₃        s₃        s₃     s₃      s₃  ⎞\n",
       "β⋅p₀⋅⎜-s₃ + ──── - ─────── + ─────── - ─────── + ─────── + ── - ───────⎟\n",
       "     ⎝      θ₂ ₀   s₂⋅θ₂ ₀   s₂⋅θ₁ ₀   s₁⋅θ₁ ₀   s₁⋅θ₀ ₀   s₀   s₀⋅θ₀ ₀⎠"
      ]
     },
     "execution_count": 10,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "beta = symbols('beta')\n",
    "p_0 = symbols('p_0')\n",
    "exec('output = p_0 * beta * expand(s_simp / s_%d)' % depth_of_tree)\n",
    "output"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 73,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAN8AAAAyBAMAAAAuDRNDAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAIol2mRBmq1TdRDK7ze+zDnIKAAAACXBIWXMAAA7EAAAOxAGVKw4bAAAEmUlEQVRYCaVY3YsbVRT/JesknUx2E+ofsKHBh4KyYbN0EWkbkD70RQf0RWjdiKCI2I0+7CKCRvTFh3b3wSqxL0HFF1sa3RdFcFeEFfGhgSrigxgrorBFti4rZf1Yz73JzNwzkzNtZy7M3HPO73xkzpx77twAwsh3BSCx+JQba7oUiyYBne/irPKVODQZ9kac2TNxYELsrpZsmPlbxpIjb8mmdkXGkiPn5LI5k9xrjOVkQwSviYgBVGflX0xq2fo9hrIiCzdCAp+d2PVJmSj2Cg0ZBb7B+2H427DA4+2YevJ0MI1iM3NkvucLOJHZxibufvIlU7oipWRd8mJaT2+Wka9MDEyZQWf+eQJ4AccMERZqJmfQLxu0SGb3/4TTy/YlhZn9Cj7DfSZcXDU5g/7NoGWycMUF7KaocFFVAlt7OaFqnL9EJwGQXUXJxcTVQMKpwwAFPMqSaAmOs8IPYR4na3idBKeYR0PjKg40AOtLQwQIqSsOmNZ4xqnXm4Q4UhNcqx9Rhrz+vlYiHJw7oWf/NrXhk7cgHu9lhCxpQ2sL003TxYxiMk08bQqBhQbnZc52c7/KKHAeZ6iwgrHYJfoRILRTrUvvRVl+oG7ecKrzZY8eN39cf56JtWfyvsykWGlznnGvMe4OGZ27Ggrnud1ij/OMSxWw1KIW/tFXP7p494u1eddzfLzsUXrmbT9dwA3aRfDUpz1UXmxniRmOGT+05nnbTxVwqq8CItfPtL9Hbtt6bhjwYRYw1PZTBbRXdZEWti3cRHZ37f5hwFf0dLqjxgWYbf9kp/NHp/O2wvf9QYz1g1buccDXUMS2wuwB8CFAW0xhF9RfRstjGFAp6MHbfronHAAVoNgCtWIq2VFAnlKAtf10ASmlfeDzMuj75pjrBZwpDx9teD/M236qgFQ01hIK71A323C2/I5zvGcGDLX9dAFbOPBJvUruV96bo8capZQv/FDbTxWw1FCrQo2T+j4KuNLVnHEz2v4dBZxl31CUyIoX8E0z4HrTiDUkg7bPmndEjwuKzce4gJo37RU0rJvqfvCnB9UEevDwCLf9MD6evxclVg5qW6BlSOPsv009D29232BSkM4eSm1mH+rSHpYfeFS6Ob+K0NsJtRTPfe52vvQB8WzhAQs1hOrvuheCz45+o1wW5cSzhQ/QeuY9JLMTdaMlwtcc11ZnCxya5ULF+cCV5eX/kK3SUqBLDelDGA9p+BY3dbbI1gobETUfuA5rB1Vc1JdSo5c6flzujpczqTpbTPUQ/Wj2gR1QNWyhVFaXsi3VmIeAsRsBLVN0tpgu4/eowgiw9mA3Mjcw9SpdbaW26EaVteR2ypQ2tJK77vLTirIOgD2cc51d2M/S1VSQ3rkVERnXIpKIQJ8tKOAvYSQAlqwLoIBTD9ClkikfueVnD7zrswWlNLKEAuBR2oEopSVKqe44xVZgHqImKyGBwI4tGqUbAEbRXJZeIZ2JfhYihMS0LPoh0ZANgCou0bK4pMRj6ss3PupT8cShOeFX+0C2eoIWvjqf5eWMEliJD5QEPRtrtBmLJgHj/75EsZvEaZzNaZ38/wHprjeknlSocAAAAABJRU5ErkJggg==\n",
      "text/latex": [
       "$$\\beta p_{1} \\left(- s_{3} + \\frac{s_{3}}{s_{0}} - \\frac{s_{3}}{s_{0} \\theta_{0 0}}\\right)$$"
      ],
      "text/plain": [
       "     ⎛      s₃      s₃  ⎞\n",
       "β⋅p₁⋅⎜-s₃ + ── - ───────⎟\n",
       "     ⎝      s₀   s₀⋅θ₀ ₀⎠"
      ]
     },
     "execution_count": 73,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#### Deriving the cross-elasticity formular for vehicles in the same make_region (second to the lowest level)\n",
    "level_of_similarity = 0 # starts from 0\n",
    "diff_node = total_end_nodes // number_splits_non_end_node ** level_of_similarity - 1\n",
    "exec('s_diff_v_x = diff(s, v_%d)' % diff_node)\n",
    "s_v_0 = s\n",
    "s_simp = s_diff_v_x.subs(sum_exp_0_0 / (1+sum_exp_0_0), s_0)\n",
    "s_v_0 = s_v_0.subs(sum_exp_0_0 / (1+sum_exp_0_0), s_0)\n",
    "for i in range(level_of_similarity):\n",
    "    exec('s_simp = s_simp.subs(sum_exp_%d_0 / sum_%d_0, s_%d / s_%d)' % (i+1, i, i+1, i))\n",
    "    exec('s_v_0 = s_v_0.subs(sum_exp_%d_0 / sum_%d_0, s_%d / s_%d)' % (i+1, i, i+1, i))\n",
    "for i in range(level_of_similarity, depth_of_tree):\n",
    "    exec('s_simp = s_simp.subs(sum_exp_%d_%d / sum_%d_%d, s_%d / s_%d)' % (i+1, diff_node // number_end_group // number_splits_non_end_node ** (depth_of_tree-i-2),\n",
    "                                                                           i, diff_node // number_end_group // number_splits_non_end_node ** (depth_of_tree-i-1), i+1, i))\n",
    "    exec('s_v_0 = s_v_0.subs(sum_exp_%d_%d / sum_%d_%d, s_%d / s_%d)' % (i+1, diff_node // number_end_group // number_splits_non_end_node ** (depth_of_tree-i-2),\n",
    "                                                                           i, diff_node // number_end_group // number_splits_non_end_node ** (depth_of_tree-i-1), i+1, i))\n",
    "exec('s_simp = s_simp.subs(exp(v_%d/theta_%d_%d) / sum_%d_%d, s_%d / s_%d)' % (diff_node, depth_of_tree-1, diff_node // number_end_group,\n",
    "                                                                           depth_of_tree-1, diff_node // number_end_group, depth_of_tree,\n",
    "                                                                              depth_of_tree-1))\n",
    "exec('s_v_0 = s_v_0.subs(exp(v_%d/theta_%d_%d) / sum_%d_%d, s_%d / s_%d)' % (diff_node, depth_of_tree-1, diff_node // number_end_group,\n",
    "                                                                           depth_of_tree-1, diff_node // number_end_group, depth_of_tree,\n",
    "                                                                              depth_of_tree-1))\n",
    "beta = Symbol('beta')\n",
    "p_1 = Symbol('p_1')\n",
    "beta * p_1 * simplify(s_simp / s_v_0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 79,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAEIAAAATBAMAAADWqMaJAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAEM3dMiKJdplmq1REu+/l2FR3AAAACXBIWXMAAA7EAAAOxAGVKw4bAAABOUlEQVQoFZWRO0sDQRDH//f2LhdyiI2CGiu18sB0NvsJJIVNEMmhxVXBVDlMIWkEy0PsjX6CFAEthdhba6FJ4TcQiaDo7MMDLRacYp6/3Z2ZBf4hTnoa6fEG/L6WMNooxVrCAwL9HUSM9H3c7V2tIGztZJfqqbTz50Din9fgBltoSaLMwp5ilUmAM+yafVxgvkq5Jsq5c9RlmN3kUgONgntElQEO5zI+U/OhjiDxxuRKCRlwTFCEJzm1+/UGm7lVVQdc8h6BAxhTtZfwmTq18oIwAfsFWIM3EYQboxLBWy0ALOQY3cL4xDYThNnGDVUXqT0lw+sOBaWP7IR0TBdmGX/Afv0BkAhvhmqSEOESc96Fw1UiPKvHTfGDVuRPRJ4UXwfJ8oCUsz+VEey0Wxd5UnwddHhjnZtf8g2Y/zne7NvK1wAAAABJRU5ErkJggg==\n",
      "text/latex": [
       "$$- \\beta p_{1} s_{3}$$"
      ],
      "text/plain": [
       "-β⋅p₁⋅s₃"
      ]
     },
     "execution_count": 79,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Cross-elasticity with outside option\n",
    "outside_s = 1 - sum_exp_0_0 / (1+sum_exp_0_0)\n",
    "s_diff_v_0 = diff(outside_s, v_0)\n",
    "s_simp = s_diff_v_0 / outside_s\n",
    "s_simp = s_simp.subs(sum_exp_0_0 / (1+sum_exp_0_0), s_0)\n",
    "for i in range(depth_of_tree):\n",
    "    exec('s_simp = s_simp.subs(sum_exp_%d_0 / sum_%d_0, s_%d / s_%d)' % (i+1, i, i+1, i))\n",
    "beta = Symbol('beta')\n",
    "p_1 = Symbol('p_1')\n",
    "beta * p_1 * simplify(s_simp)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
