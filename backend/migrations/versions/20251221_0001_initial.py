from alembic import op
import sqlalchemy as sa
import uuid
from sqlalchemy.dialects import postgresql

revision = "20251221_0001"
down_revision = None

def upgrade():
    op.create_table(
        "molecules",
        sa.Column("id", postgresql.UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
        sa.Column("smiles", sa.String(), nullable=False),
        sa.Column("name", sa.String(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.text("NOW()")),
    )
    op.create_index("ix_molecules_smiles", "molecules", ["smiles"])

def downgrade():
    op.drop_index("ix_molecules_smiles", table_name="molecules")
    op.drop_table("molecules")
